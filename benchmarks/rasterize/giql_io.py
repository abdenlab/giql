# giql_io.py
# Author: Conrad Bzura <conradbzura@gmail.com>

"""GIQL / polars-bio pre-processing engine.

Vendored into the benchmark so it runs against **unmodified upstream
cherimoya**. It imports nothing from cherimoya, only tangermeme, polars-bio
and giql, so it stands alone; `peak_generator` in `pileup_bench.py` is what
lets an alternative extractor reach the sampler without cherimoya needing an
injection point of its own.

This module is a drop-in replacement for the ``bam2bw`` + ``tangermeme.io.
extract_loci`` half of the training data path. Alignments and sequences are
read through the datafusion-bio table providers bundled with polars-bio, the
Tn5 cut-site pileup is a DataFusion aggregate over the BAM table, and the
per-locus signal and exclusion-list tests are GIQL ``INTERSECTS`` joins that
polars-bio plans as ``IntervalJoinExec``. No bigWig is ever written: the
pileup is materialized once as a sorted ``{name}.cuts.parquet`` table of
``(chrom, pos, value)`` rows and read back per chromosome.

Everything here reproduces the reference semantics exactly (see
``tests/test_giql_io.py`` and ``scripts/verify_ingestion.py``):

- every alignment record except unmapped ones contributes one cut site
  (forward: ``start + pos_shift``; reverse: ``end - 1 + neg_shift``), reads on
  contigs absent from the FASTA are dropped, and all input files collapse into
  one channel, exactly like a single ``bam2bw`` invocation;
- loci are loaded, filtered by chromosome and interleaved with tangermeme's
  own loader, windows use tangermeme's centre / bounds / jitter arithmetic,
  and the exclusion test reproduces its 100-bp bin masks;
- sequences are one-hot encoded with ``tangermeme.utils.one_hot_encode``.

polars-bio holds one process-global session; all tables registered here use a
``cherimoya_`` prefix so they never collide with a caller's tables.
"""

import collections
import os

import numpy
import pandas
import torch

CUT_SITE_SUFFIX = ".cuts.parquet"

# polars-bio's BAM provider returns ``start`` equal to pysam's 0-based
# ``reference_start`` and ``end`` equal to the exclusive ``reference_end`` on
# the default session, despite its docstring claiming 1-based coordinates.
# ``tests/test_giql_io.py`` pins this against pysam; bump the offset (and the
# test) if a future polars-bio release changes the convention.
_BAM_START_OFFSET = 0

_FASTA_TABLE = "cherimoya_fasta"
_SIZES_TABLE = "cherimoya_chrom_sizes"
_WINDOWS_TABLE = "cherimoya_windows"
_CUTS_TABLE = "cherimoya_cuts"
_LOCUS_BINS_TABLE = "cherimoya_locus_bins"
_EXCLUSION_BINS_TABLE = "cherimoya_exclusion_bins"

_cache = {
    "fasta": None,
    "chrom_lengths": None,
    "sequences": {},
    "sql": {},
}


def clear_cache():
    """Forget cached chromosome lengths, sequences and transpiled queries."""

    _cache["fasta"] = None
    _cache["chrom_lengths"] = None
    _cache["sequences"] = {}
    _cache["sql"] = {}


###
# Engine access
###


def _pb():
    """Import polars-bio lazily with quiet defaults."""

    os.environ.setdefault("TQDM_DISABLE", "1")
    import polars_bio as pb

    try:
        pb.set_loglevel("error")
    except Exception:
        pass
    return pb


def _giql_sql(key, query, tables):
    """Transpile a GIQL query for the datafusion-bio target once per process."""

    if key not in _cache["sql"]:
        from giql import transpile

        _cache["sql"][key] = transpile(query, tables=tables, dialect="datafusion-bio")
    return _cache["sql"][key]


def _interval_tables(*specs):
    from giql import Table

    return [
        Table(name, chrom_col="chrom", start_col=start, end_col=end, strand_col=None)
        for name, start, end in specs
    ]


def _register_frame(pb, name, frame):
    pb.from_polars(name, frame)


def _explain(pb, sql):
    """Return the physical plan of ``sql`` as one string."""

    import polars

    plan = pb.sql("EXPLAIN " + sql).collect()
    return (
        "\n".join(
            str(v)
            for col in plan.get_columns()
            for v in col.to_list()
            if isinstance(v, str)
        )
        if isinstance(plan, polars.DataFrame)
        else str(plan)
    )


###
# Sequences
###


def _register_fasta(sequences):
    """Register ``sequences`` (a FASTA path) once, dropping stale caches."""

    if not isinstance(sequences, str):
        raise NotImplementedError(
            "The GIQL engine reads sequences from a "
            "FASTA file; in-memory sequence dictionaries are not supported."
        )

    pb = _pb()
    if _cache["fasta"] != sequences:
        pb.register_fasta(sequences, _FASTA_TABLE)
        _cache["fasta"] = sequences
        _cache["chrom_lengths"] = None
        _cache["sequences"] = {}
    return pb


def chrom_lengths(sequences):
    """Return an ordered ``chrom -> length`` mapping for a FASTA file.

    The lengths come from one scan of the polars-bio FASTA table; pyfaidx is
    the fallback (and the reference definition used by bam2bw and tangermeme)
    if the provider cannot hold a genome in a single query.
    """

    pb = _register_fasta(sequences)
    if _cache["chrom_lengths"] is None:
        try:
            frame = pb.sql(
                "SELECT name, LENGTH(sequence) AS length FROM {}".format(_FASTA_TABLE)
            ).collect()
            lengths = collections.OrderedDict(
                zip(frame["name"].to_list(), frame["length"].to_list())
            )
        except Exception:
            import pyfaidx

            fasta = pyfaidx.Fasta(sequences)
            lengths = collections.OrderedDict(
                (chrom, len(fasta[chrom])) for chrom in fasta.keys()
            )
            fasta.close()

        _cache["chrom_lengths"] = lengths
    return _cache["chrom_lengths"]


def _prefetch_sequences(sequences, chroms):
    """Load every chromosome in ``chroms`` with a single FASTA scan."""

    pb = _register_fasta(sequences)
    missing = [str(c) for c in chroms if str(c) not in _cache["sequences"]]
    if not missing:
        return

    frame = pb.sql(
        "SELECT name, sequence FROM {} WHERE name IN ({})".format(
            _FASTA_TABLE, ", ".join("'{}'".format(c) for c in missing)
        )
    ).collect()
    for name, sequence in zip(frame["name"].to_list(), frame["sequence"].to_list()):
        _cache["sequences"][name] = sequence


def _chrom_sequence(sequences, chrom):
    """Return the residues of one chromosome as a Python string (cached)."""

    chrom = str(chrom)
    if chrom not in _cache["sequences"]:
        _prefetch_sequences(sequences, [chrom])
    if chrom not in _cache["sequences"]:
        raise KeyError(chrom)
    return _cache["sequences"][chrom]


###
# Cut-site pileup
###


def register_alignments(paths, prefix="cherimoya_bam"):
    """Register BAM/SAM files as polars-bio tables and return their names."""

    pb = _pb()
    names = []
    for i, path in enumerate(paths):
        if not path.endswith((".bam", ".sam")):
            raise ValueError(
                "The GIQL engine converts BAM/SAM alignments; got {}".format(path)
            )

        name = "{}_{}".format(prefix, i)
        if path.endswith(".sam"):
            pb.register_sam(path, name)
        else:
            pb.register_bam(path, name)
        names.append(name)

    return names


def pileup_sql(
    bam_tables,
    pos_shift=0,
    neg_shift=0,
    unstranded=False,
    scale_factor=1.0,
    sizes_table=_SIZES_TABLE,
):
    """Build the cut-site aggregate over one or more registered BAM tables.

    Mirrors bam2bw: forward reads (``flags & 16 = 0``) cut at ``start +
    pos_shift``; reverse reads cut at ``end - 1 + neg_shift``; unmapped records
    (``flags & 4``) are skipped; records on contigs absent from the sizes
    table are dropped; several files are concatenated before counting. The
    one deliberate difference is that cut sites shifted outside ``[0, len)``
    are dropped here, whereas bam2bw hands them to pyBigWig unclipped; no
    training window can read those positions, so tensors are unaffected.
    """

    strand = (
        "CASE WHEN (CAST(b.flags AS BIGINT) & 16) = 0 THEN '+' ELSE '-' END AS strand"
    )
    pos = (
        "CASE WHEN (CAST(b.flags AS BIGINT) & 16) = 0 "
        "THEN CAST(b.start AS BIGINT) + {offset} + ({pos_shift}) "
        'ELSE CAST(b."end" AS BIGINT) + {offset} + ({neg_shift}) - 1 END AS pos'
    ).format(
        offset=_BAM_START_OFFSET, pos_shift=int(pos_shift), neg_shift=int(neg_shift)
    )

    branches = [
        "SELECT b.chrom, {pos}, {strand} FROM {table} AS b "
        "WHERE (CAST(b.flags AS BIGINT) & 4) = 0".format(
            pos=pos, strand=strand, table=table
        )
        for table in bam_tables
    ]

    group = "c.chrom, c.pos" if unstranded else "c.chrom, c.pos, c.strand"
    return (
        "WITH cuts AS ({branches}) "
        "SELECT {group}, CAST(COUNT(*) AS DOUBLE) * {scale} AS value "
        "FROM cuts AS c JOIN {sizes} AS s ON s.chrom = c.chrom "
        "WHERE c.pos >= 0 AND c.pos < s.length "
        "GROUP BY {group}"
    ).format(
        branches=" UNION ALL ".join(branches),
        group=group,
        scale=float(scale_factor),
        sizes=sizes_table,
    )


def _cut_site_paths(name, unstranded):
    if unstranded:
        return [name + CUT_SITE_SUFFIX]
    return [name + ".+" + CUT_SITE_SUFFIX, name + ".-" + CUT_SITE_SUFFIX]


def pileup_alignments(
    alignments,
    sequences,
    name,
    pos_shift=0,
    neg_shift=0,
    unstranded=False,
    fragments=False,
    scale_factor=1.0,
    read_depth=False,
    verbose=False,
):
    """Convert BAM/SAM files into ``{name}.cuts.parquet`` cut-site tables.

    Parameters
    ----------
    alignments: list of str
            BAM/SAM paths. All files are pooled into the same channel(s).

    sequences: str
            FASTA path. Its contigs and lengths define which reads are kept.

    name: str
            Output prefix. Unstranded data yields ``{name}.cuts.parquet``;
            stranded data yields ``{name}.+.cuts.parquet`` and
            ``{name}.-.cuts.parquet``.

    pos_shift, neg_shift: int
            Shifts applied to forward-strand 5' ends and to reverse-strand 5'
            ends, respectively (bam2bw's ``-ps`` / ``-ns``).

    unstranded: bool
            Pool both strands into a single table.

    fragments, read_depth: bool
            Not supported (the ATAC-seq BAM recipe never sets them).

    scale_factor: float
            Multiplier applied to every count.

    Returns
    -------
    paths: list of str
            The written Parquet paths, in channel order.
    """

    if fragments:
        raise NotImplementedError(
            "fragments=True is not supported by the "
            "GIQL engine; use the bam2bw engine for fragment files."
        )
    if read_depth:
        raise NotImplementedError("read_depth=True is not supported by the GIQL engine.")

    import polars

    pb = _pb()
    lengths = chrom_lengths(sequences)
    tables = register_alignments(alignments)

    _register_frame(
        pb,
        _SIZES_TABLE,
        polars.DataFrame(
            {
                "chrom": list(lengths.keys()),
                "length": list(lengths.values()),
            },
            schema={"chrom": polars.Utf8, "length": polars.Int64},
        ),
    )

    try:
        pb.set_option(
            "datafusion.execution.target_partitions", str(max(1, os.cpu_count() or 1))
        )
    except Exception:
        pass

    sql = pileup_sql(
        tables,
        pos_shift=pos_shift,
        neg_shift=neg_shift,
        unstranded=unstranded,
        scale_factor=scale_factor,
    )
    if verbose:
        print("GIQL engine: piling up cut sites from {}".format(", ".join(alignments)))

    frame = pb.sql(sql).collect()
    frame = frame.with_columns(
        polars.col("pos").cast(polars.Int64),
        polars.col("value").cast(polars.Float64),
    )

    paths = _cut_site_paths(name, unstranded)
    if unstranded:
        parts = [frame]
    else:
        parts = [
            frame.filter(polars.col("strand") == s).drop("strand") for s in ("+", "-")
        ]

    for path, part in zip(paths, parts):
        part = part.select("chrom", "pos", "value").sort("chrom", "pos")
        part.write_parquet(path, row_group_size=1_000_000, statistics=True)
        if verbose:
            print("GIQL engine: wrote {} cut sites to {}".format(part.height, path))

    return paths


###
# Loci
###


def load_windows(
    loci,
    chrom_lengths,
    chroms=None,
    summits=False,
    in_window=2114,
    out_window=1000,
    max_jitter=0,
    has_signals=True,
):
    """Load loci the way tangermeme does and derive every window bound.

    Returns a DataFrame with one row per locus (after chromosome filtering
    and interleaving, in tangermeme's order) and columns ``chrom``, ``mid``,
    ``ws``/``we`` (the bounds window), ``in_bounds``, ``sig_start``/``sig_end``
    and ``seq_start``/``seq_end``.
    """

    from tangermeme.io import _interleave_loci

    loci = _interleave_loci(loci, chroms, summits=summits)

    in_width, out_width = in_window // 2, out_window // 2
    if not has_signals:
        out_width = 0

    max_width = max(in_width, out_width)

    start = loci["start"].to_numpy().astype(numpy.int64)
    end = loci["end"].to_numpy().astype(numpy.int64)
    mid = start + (end - start) // 2

    lengths = loci["chrom"].map(chrom_lengths)
    if lengths.isnull().any():
        raise KeyError(loci["chrom"][lengths.isnull()].iloc[0])
    lengths = lengths.to_numpy().astype(numpy.int64)

    ws = mid - max_width - max_jitter
    we = mid + max_width + max_jitter

    return pandas.DataFrame(
        {
            "chrom": loci["chrom"].to_numpy(),
            "mid": mid,
            "ws": ws,
            "we": we,
            # tangermeme >= 1.4 keeps a window whose end equals the chromosome
            # length (``end > chrom_len`` drops); older releases used ``>=``.
            "in_bounds": (ws >= 0) & (we <= lengths),
            "sig_start": mid - out_width - max_jitter,
            "sig_end": mid + out_width + max_jitter + (out_window % 2),
            "seq_start": mid - in_width - max_jitter,
            "seq_end": mid + in_width + max_jitter + (in_window % 2),
        }
    )


def exclusion_mask(windows, exclusion_lists):
    """Flag loci whose bounds window touches an excluded 100-bp bin.

    tangermeme marks bins ``start // 100`` through ``end // 100`` (inclusive)
    for every exclusion row and drops a locus when any bin from ``ws // 100``
    through ``we // 100`` is marked. Both are half-open bin intervals, so
    the test is a GIQL ``INTERSECTS`` join between the two bin tables.
    """

    import polars

    names = ("chrom", "start", "end")
    exclusions = pandas.concat(
        [
            pandas.read_csv(path, sep="\t", names=names, header=None, usecols=(0, 1, 2))
            for path in exclusion_lists
        ]
    )

    mask = numpy.zeros(len(windows), dtype=bool)
    candidates = numpy.flatnonzero(windows["in_bounds"].to_numpy())
    if len(exclusions) == 0 or len(candidates) == 0:
        return mask

    pb = _pb()
    _register_frame(
        pb,
        _EXCLUSION_BINS_TABLE,
        polars.DataFrame(
            {
                "chrom": exclusions["chrom"].astype(str).to_numpy(),
                "bstart": exclusions["start"].to_numpy().astype(numpy.int64) // 100,
                "bend": exclusions["end"].to_numpy().astype(numpy.int64) // 100 + 1,
            },
            schema={"chrom": polars.Utf8, "bstart": polars.Int64, "bend": polars.Int64},
        ),
    )

    _register_frame(
        pb,
        _LOCUS_BINS_TABLE,
        polars.DataFrame(
            {
                "idx": candidates.astype(numpy.int64),
                "chrom": windows["chrom"].to_numpy()[candidates].astype(str),
                "bstart": windows["ws"].to_numpy()[candidates] // 100,
                "bend": windows["we"].to_numpy()[candidates] // 100 + 1,
            },
            schema={
                "chrom": polars.Utf8,
                "idx": polars.Int64,
                "bstart": polars.Int64,
                "bend": polars.Int64,
            },
        ),
    )

    sql = _giql_sql(
        "exclusion",
        "SELECT DISTINCT l.idx FROM {} AS l JOIN {} AS e "
        "ON l.interval INTERSECTS e.interval".format(
            _LOCUS_BINS_TABLE, _EXCLUSION_BINS_TABLE
        ),
        _interval_tables(
            (_LOCUS_BINS_TABLE, "bstart", "bend"),
            (_EXCLUSION_BINS_TABLE, "bstart", "bend"),
        ),
    )

    excluded = pb.sql(sql).collect()["idx"].to_numpy()
    mask[excluded] = True
    return mask


def _cut_site_frame(source):
    """Normalize one cut-site source to a polars LazyFrame.

    A source is either a path to a ``.cuts.parquet`` table written by
    :func:`pileup_alignments`, or an in-memory polars ``DataFrame`` /
    ``LazyFrame`` with the same ``chrom, pos, value`` columns. The in-memory
    form is what lets a caller run BAM to tensors without touching disk.
    """

    import polars

    if isinstance(source, str):
        return polars.scan_parquet(source)
    if isinstance(source, polars.LazyFrame):
        return source
    if isinstance(source, polars.DataFrame):
        return source.lazy()

    raise NotImplementedError(
        "The GIQL engine reads signals from {} tables "
        "or in-memory polars frames; got {!r}".format(
            CUT_SITE_SUFFIX, type(source).__name__
        )
    )


def signal_sql():
    """The transpiled per-locus signal join (exposed for plan inspection)."""

    return _giql_sql(
        "signal",
        "SELECT w.idx, c.cstart - w.wstart AS offset, c.value "
        "FROM {} AS w JOIN {} AS c ON w.interval INTERSECTS c.interval".format(
            _WINDOWS_TABLE, _CUTS_TABLE
        ),
        _interval_tables(
            (_WINDOWS_TABLE, "wstart", "wend"), (_CUTS_TABLE, "cstart", "cend")
        ),
    )


def extract_signal(windows, keep, cut_sites, verbose=False):
    """Extract dense per-locus signal from cut-site tables.

    Parameters
    ----------
    windows: pandas.DataFrame
            As returned by :func:`load_windows`.

    keep: numpy.ndarray, dtype=bool
            Which loci to extract; output rows follow the kept loci in order.

    cut_sites: list of str or polars frames
            One channel per entry, each either a path to a ``.cuts.parquet``
            table or an in-memory polars ``DataFrame`` / ``LazyFrame`` of
            ``chrom, pos, value``.

    Returns
    -------
    signal: numpy.ndarray, shape=(keep.sum(), len(cut_sites), L), float32
            Cut-site counts per position, zero where no cut site exists.
    """

    import polars

    kept = windows[keep].reset_index(drop=True)
    length = int((kept["sig_end"] - kept["sig_start"]).iloc[0]) if len(kept) else 0
    signal = numpy.zeros((len(kept), len(cut_sites), length), dtype=numpy.float32)
    if len(kept) == 0:
        return signal

    pb = _pb()
    sql = signal_sql()
    channels = [_cut_site_frame(source) for source in cut_sites]
    chroms = pandas.unique(kept["chrom"])
    for chrom in chroms:
        rows = numpy.flatnonzero((kept["chrom"] == chrom).to_numpy())
        _register_frame(
            pb,
            _WINDOWS_TABLE,
            polars.DataFrame(
                {
                    "idx": rows.astype(numpy.int64),
                    "chrom": numpy.repeat(str(chrom), len(rows)),
                    "wstart": kept["sig_start"].to_numpy()[rows].astype(numpy.int64),
                    "wend": kept["sig_end"].to_numpy()[rows].astype(numpy.int64),
                },
                schema={
                    "idx": polars.Int64,
                    "chrom": polars.Utf8,
                    "wstart": polars.Int64,
                    "wend": polars.Int64,
                },
            ),
        )

        for channel, frame in enumerate(channels):
            cuts = (
                frame.filter(polars.col("chrom") == str(chrom))
                .select(
                    polars.col("chrom").cast(polars.Utf8),
                    polars.col("pos").cast(polars.Int64).alias("cstart"),
                    (polars.col("pos").cast(polars.Int64) + 1).alias("cend"),
                    polars.col("value").cast(polars.Float64),
                )
                .collect()
            )
            if cuts.height == 0:
                continue

            _register_frame(pb, _CUTS_TABLE, cuts)
            hits = pb.sql(sql).collect()
            if hits.height:
                signal[hits["idx"].to_numpy(), channel, hits["offset"].to_numpy()] = (
                    hits["value"].to_numpy().astype(numpy.float32)
                )

        if verbose:
            print(
                "GIQL engine: extracted signal for {} loci on {}".format(
                    len(rows), chrom
                )
            )

    return signal


def extract_sequences(
    windows, keep, sequences, alphabet=["A", "C", "G", "T"], ignore=["N"], verbose=False
):
    """One-hot encode the sequence window of every kept locus."""

    from tangermeme.utils import one_hot_encode

    kept = windows[keep]
    _prefetch_sequences(sequences, pandas.unique(kept["chrom"]))

    seqs = []
    current, seq = None, None
    for chrom, start, end in zip(kept["chrom"], kept["seq_start"], kept["seq_end"]):
        if chrom != current:
            current, seq = chrom, _chrom_sequence(sequences, chrom)

        seqs.append(
            one_hot_encode(
                seq[int(start) : int(end)].upper(), alphabet=alphabet, ignore=ignore
            )
        )

    if len(seqs) == 0:
        return torch.zeros((0, len(alphabet), 0), dtype=torch.int8)
    return torch.from_numpy(numpy.stack(seqs))


def extract_loci(
    loci,
    sequences,
    signals=None,
    in_signals=None,
    chroms=None,
    in_window=2114,
    out_window=1000,
    max_jitter=0,
    min_counts=None,
    max_counts=None,
    target_idx=0,
    n_loci=None,
    summits=False,
    alphabet=["A", "C", "G", "T"],
    ignore=["N"],
    exclusion_lists=None,
    return_mask=False,
    verbose=False,
):
    """Extract sequence and signal for each locus (tangermeme-compatible).

    This has the signature, return structure, ordering and dtypes of
    ``tangermeme.io.extract_loci`` so it can be handed to
    :func:`cherimoya.io.PeakGenerator` as ``extractor``. Differences:
    ``sequences`` must be a FASTA path, each entry of ``signals`` must be a
    ``.cuts.parquet`` path produced by :func:`pileup_alignments` or an
    equivalent in-memory polars frame, and ``in_signals`` (controls) and
    ``n_loci`` are not supported.
    """

    if in_signals is not None:
        raise NotImplementedError(
            "Control tracks (in_signals) are not supported by the GIQL engine."
        )
    if n_loci is not None:
        raise NotImplementedError("n_loci is not supported by the GIQL engine.")
    if signals is not None:
        if not isinstance(signals, (list, tuple)):
            raise ValueError("Provided signals must be in the form of a list.")
        for signal in signals:
            if isinstance(signal, str) and not signal.endswith(CUT_SITE_SUFFIX):
                raise NotImplementedError(
                    "The GIQL engine reads signals from {} tables; got {!r}".format(
                        CUT_SITE_SUFFIX, signal
                    )
                )
            if not isinstance(signal, str):
                # Raises for anything that is not a polars frame.
                _cut_site_frame(signal)

    lengths = chrom_lengths(sequences)
    windows = load_windows(
        loci,
        lengths,
        chroms=chroms,
        summits=summits,
        in_window=in_window,
        out_window=out_window,
        max_jitter=max_jitter,
        has_signals=signals is not None,
    )

    keep = windows["in_bounds"].to_numpy().copy()
    if exclusion_lists is not None:
        keep &= ~exclusion_mask(windows, exclusion_lists)

    signal = None
    if signals is not None:
        signal = extract_signal(windows, keep, signals, verbose=verbose)

        if min_counts is not None or max_counts is not None:
            totals = signal[:, target_idx].sum(axis=1)
            ok = numpy.ones(len(totals), dtype=bool)
            if min_counts is not None:
                ok &= ~(totals < min_counts)
            if max_counts is not None:
                ok &= ~(totals > max_counts)

            keep[numpy.flatnonzero(keep)[~ok]] = False
            signal = signal[ok]

    seqs = extract_sequences(
        windows, keep, sequences, alphabet=alphabet, ignore=ignore, verbose=verbose
    )

    y_return = [seqs]
    if signals is not None:
        y_return.append(torch.from_numpy(signal))
    if return_mask:
        y_return.append(torch.from_numpy(keep))

    return y_return[0] if len(y_return) == 1 else y_return
