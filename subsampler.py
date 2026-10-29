#!/usr/bin/env python3

import argparse
import gzip
import random
import subprocess
import sys
from pathlib import Path

import requests


# ============================================================
# ENA CONFIGURATION
# ============================================================

ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"


# ============================================================
# ENA FUNCTIONS
# ============================================================

def get_runs(accession):
    """
    Return run accessions for an ENA/SRA/DRA accession.

    If accession is already a run accession (SRR/ERR/DRR),
    return it directly.
    """

    prefixes = ("SRR", "ERR", "DRR")

    if accession.startswith(prefixes):
        return [accession]

    params = {
        "accession": accession,
        "result": "read_run",
        "fields": "run_accession",
        "format": "tsv",
    }

    response = requests.get(
        ENA_API,
        params=params,
        timeout=60,
    )

    response.raise_for_status()

    lines = response.text.strip().splitlines()

    if len(lines) < 2:
        return []

    return [
        line.split("\t")[0]
        for line in lines[1:]
        if line.strip()
    ]


def fetch_fastq(accession):
    """
    Retrieve FASTQ FTP URLs for a run accession.
    """

    params = {
        "accession": accession,
        "result": "read_run",
        "fields": "fastq_ftp",
        "format": "tsv",
    }

    response = requests.get(
        ENA_API,
        params=params,
        timeout=60,
    )

    response.raise_for_status()

    lines = response.text.strip().splitlines()

    if len(lines) < 2:
        return []

    fields = lines[1].split("\t")

    if len(fields) < 2:
        return []

    urls = fields[1]

    if not urls:
        return []

    return [
        "https://" + url if not url.startswith("http") else url
        for url in urls.split(";")
        if url
    ]


def download_file(url, output_file):
    """
    Download a single FASTQ file.

    Downloads to a temporary file first and only renames it
    after a successful download.
    """

    output_file = Path(output_file)
    temporary_file = Path(str(output_file) + ".download")

    if temporary_file.exists():
        temporary_file.unlink()

    print(f"  Downloading:")
    print(f"    {url}")
    print(f"  -> {output_file}")

    try:
        with open(temporary_file, "wb") as out:
            subprocess.run(
                [
                    "curl",
                    "-L",
                    "--fail",
                    "--retry",
                    "3",
                    "--retry-delay",
                    "5",
                    url,
                ],
                stdout=out,
                check=True,
            )

    except subprocess.CalledProcessError:
        if temporary_file.exists():
            temporary_file.unlink()

        raise

    if (
        not temporary_file.exists()
        or temporary_file.stat().st_size == 0
    ):
        if temporary_file.exists():
            temporary_file.unlink()

        raise RuntimeError(
            f"Downloaded file is empty: {output_file}"
        )

    temporary_file.rename(output_file)


def download_sample(accession, raw_reads_dir):
    """
    Download all FASTQ files belonging to an accession.

    For paired-end data:
        SAMPLE_R1.fastq.gz
        SAMPLE_R2.fastq.gz

    For single-end data:
        SAMPLE.fastq.gz

    Returns a list of downloaded run accessions.
    """

    raw_reads_dir = Path(raw_reads_dir)
    raw_reads_dir.mkdir(parents=True, exist_ok=True)

    runs = get_runs(accession)

    if not runs:
        print(f"[DOWNLOAD FAILED] No runs found for {accession}")
        return []

    downloaded_runs = []

    for run in runs:

        # ----------------------------------------------------
        # Check whether this run already exists locally
        # ----------------------------------------------------

        existing_r1 = raw_reads_dir / f"{run}_R1.fastq.gz"
        existing_r2 = raw_reads_dir / f"{run}_R2.fastq.gz"
        existing_single = raw_reads_dir / f"{run}.fastq.gz"

        if existing_r1.exists() and existing_r2.exists():
            print(f"[LOCAL] {run}: paired FASTQ already exists")
            downloaded_runs.append(run)
            continue

        if existing_single.exists():
            print(f"[LOCAL] {run}: FASTQ already exists")
            downloaded_runs.append(run)
            continue

        # ----------------------------------------------------
        # Ask ENA for FASTQ files
        # ----------------------------------------------------

        try:
            urls = fetch_fastq(run)

        except requests.exceptions.RequestException as exc:
            print(f"[DOWNLOAD FAILED] ENA query failed for {run}")
            print(f"  {exc}")
            continue

        if not urls:
            print(f"[DOWNLOAD FAILED] No FASTQ files found for {run}")
            continue

        try:

            # ------------------------------------------------
            # Paired-end
            # ------------------------------------------------

            if len(urls) >= 2:

                if len(urls) > 2:
                    print(
                        f"[WARNING] {run} has {len(urls)} FASTQ URLs. "
                        f"Using the first two."
                    )

                r1 = raw_reads_dir / f"{run}_R1.fastq.gz"
                r2 = raw_reads_dir / f"{run}_R2.fastq.gz"

                download_file(urls[0], r1)
                download_file(urls[1], r2)

                downloaded_runs.append(run)

            # ------------------------------------------------
            # Single-end
            # ------------------------------------------------

            else:

                outfile = raw_reads_dir / f"{run}.fastq.gz"

                download_file(urls[0], outfile)

                downloaded_runs.append(run)

        except (
            subprocess.CalledProcessError,
            OSError,
            RuntimeError,
        ) as exc:

            print(f"[DOWNLOAD FAILED] {run}")
            print(f"  {exc}")

            # Remove incomplete files
            for file in [
                raw_reads_dir / f"{run}_R1.fastq.gz",
                raw_reads_dir / f"{run}_R2.fastq.gz",
                raw_reads_dir / f"{run}.fastq.gz",
            ]:
                if file.exists():
                    file.unlink()

    return downloaded_runs


# ============================================================
# FASTQ HELPERS
# ============================================================

def count_fastq_reads(filename):
    """
    Count reads in a FASTQ file.

    Used only for validation.
    """

    opener = gzip.open if str(filename).endswith(".gz") else open

    count = 0

    with opener(filename, "rt") as handle:
        for _ in handle:
            count += 1

    if count % 4 != 0:
        raise ValueError(
            f"Invalid FASTQ file: {filename}"
        )

    return count // 4


def validate_fastq(filename):
    """
    Check that a FASTQ file exists, is non-empty and contains
    a valid number of FASTQ lines.
    """

    filename = Path(filename)

    if not filename.exists():
        return False

    if filename.stat().st_size == 0:
        return False

    try:
        count_fastq_reads(filename)
        return True

    except Exception:
        return False


# ============================================================
# SAMPLE IDENTIFICATION
# ============================================================

def identify_samples(indir):
    """
    Identify available samples in the input directory.

    Supported paired formats:
        SAMPLE_R1.fastq.gz / SAMPLE_R2.fastq.gz
        SAMPLE_1.fastq.gz  / SAMPLE_2.fastq.gz

    Supported single-end formats:
        SAMPLE.fastq.gz
        SAMPLE.fq.gz
    """

    indir = Path(indir)

    samples = {}

    # --------------------------------------------------------
    # Paired R1/R2
    # --------------------------------------------------------

    for r1 in indir.glob("*_R1.fastq.gz"):

        sample = r1.name[:-len("_R1.fastq.gz")]
        r2 = indir / f"{sample}_R2.fastq.gz"

        if r2.exists():
            samples[sample] = {
                "type": "paired",
                "r1": r1,
                "r2": r2,
            }

    # --------------------------------------------------------
    # Paired _1/_2
    # --------------------------------------------------------

    for r1 in indir.glob("*_1.fastq.gz"):

        sample = r1.name[:-len("_1.fastq.gz")]
        r2 = indir / f"{sample}_2.fastq.gz"

        if r2.exists() and sample not in samples:
            samples[sample] = {
                "type": "paired",
                "r1": r1,
                "r2": r2,
            }

    # --------------------------------------------------------
    # Single-end
    # --------------------------------------------------------

    for fq in indir.glob("*.fastq.gz"):

        name = fq.name

        if name.endswith("_R1.fastq.gz"):
            continue

        if name.endswith("_R2.fastq.gz"):
            continue

        if name.endswith("_1.fastq.gz"):
            continue

        if name.endswith("_2.fastq.gz"):
            continue

        sample = name[:-len(".fastq.gz")]

        if sample not in samples:
            samples[sample] = {
                "type": "single",
                "r1": fq,
                "r2": None,
            }

    # --------------------------------------------------------
    # .fq.gz
    # --------------------------------------------------------

    for fq in indir.glob("*.fq.gz"):

        sample = fq.name[:-len(".fq.gz")]

        if sample not in samples:
            samples[sample] = {
                "type": "single",
                "r1": fq,
                "r2": None,
            }

    return samples


# ============================================================
# PARQUET CHECK
# ============================================================

def parquet_exists(parquet_dir, sample, reads):
    """
    Check whether the requested sample/read-depth combination
    has already been converted to Parquet.
    """

    parquet_dir = Path(parquet_dir)

    possible_files = [
        parquet_dir / f"{reads}_reads" / f"{sample}.parquet",
        parquet_dir / f"{sample}.parquet",
    ]

    return any(file.exists() for file in possible_files)


# ============================================================
# RESERVOIR SAMPLING
# ============================================================

def reservoir_sample_single(
    infile,
    outfile,
    n_reads,
    seed,
):
    """
    Reservoir-sample a single-end FASTQ file.
    """

    rng = random.Random(seed)

    reservoir = []

    opener = gzip.open if str(infile).endswith(".gz") else open

    with opener(infile, "rt") as handle:

        index = 0

        while True:

            record = [
                handle.readline(),
                handle.readline(),
                handle.readline(),
                handle.readline(),
            ]

            if not record[0]:
                break

            if not all(record):
                raise ValueError(
                    f"Incomplete FASTQ record in {infile}"
                )

            if index < n_reads:

                reservoir.append(record)

            else:

                replacement = rng.randint(0, index)

                if replacement < n_reads:
                    reservoir[replacement] = record

            index += 1

    opener_out = gzip.open if str(outfile).endswith(".gz") else open

    with opener_out(outfile, "wt") as out:

        for record in reservoir:
            out.writelines(record)

    return len(reservoir)


def reservoir_sample_paired(
    r1,
    r2,
    out1,
    out2,
    n_reads,
    seed,
):
    """
    Reservoir-sample paired-end FASTQ files while keeping
    R1/R2 pairs together.
    """

    rng = random.Random(seed)

    reservoir = []

    opener1 = gzip.open if str(r1).endswith(".gz") else open
    opener2 = gzip.open if str(r2).endswith(".gz") else open

    with opener1(r1, "rt") as h1, opener2(r2, "rt") as h2:

        index = 0

        while True:

            record1 = [
                h1.readline(),
                h1.readline(),
                h1.readline(),
                h1.readline(),
            ]

            record2 = [
                h2.readline(),
                h2.readline(),
                h2.readline(),
                h2.readline(),
            ]

            if not record1[0] and not record2[0]:
                break

            if not all(record1) or not all(record2):
                raise ValueError(
                    f"Incomplete FASTQ record in {r1} or {r2}"
                )

            if index < n_reads:

                reservoir.append(
                    (record1, record2)
                )

            else:

                replacement = rng.randint(0, index)

                if replacement < n_reads:
                    reservoir[replacement] = (
                        record1,
                        record2,
                    )

            index += 1

    opener_out1 = (
        gzip.open
        if str(out1).endswith(".gz")
        else open
    )

    opener_out2 = (
        gzip.open
        if str(out2).endswith(".gz")
        else open
    )

    with (
        opener_out1(out1, "wt") as h1,
        opener_out2(out2, "wt") as h2,
    ):

        for record1, record2 in reservoir:
            h1.writelines(record1)
            h2.writelines(record2)

    return len(reservoir)


# ============================================================
# SAMPLE PROCESSING
# ============================================================

def process_sample(
    sample,
    sample_info,
    outdir,
    reads,
    seed,
):
    """
    Perform reservoir sampling for one sample.
    """

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if sample_info["type"] == "paired":

        outfile_r1 = outdir / f"{sample}_R1.fastq.gz"
        outfile_r2 = outdir / f"{sample}_R2.fastq.gz"

        print(
            f"[SUBSAMPLING] {sample} "
            f"(paired-end, maximum {reads:,} reads)"
        )

        count = reservoir_sample_paired(
            sample_info["r1"],
            sample_info["r2"],
            outfile_r1,
            outfile_r2,
            reads,
            seed,
        )

        print(
            f"  Sampled {count:,} read pairs"
        )

        return count

    else:

        outfile = outdir / f"{sample}.fastq.gz"

        print(
            f"[SUBSAMPLING] {sample} "
            f"(single-end, maximum {reads:,} reads)"
        )

        count = reservoir_sample_single(
            sample_info["r1"],
            outfile,
            reads,
            seed,
        )

        print(
            f"  Sampled {count:,} reads"
        )

        return count


# ============================================================
# MAIN
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Download missing ENA FASTQ files and "
            "reservoir-subsample sequencing reads."
        )
    )

    parser.add_argument(
        "--indir",
        required=True,
        help="Directory containing raw FASTQ files.",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory for subsampled FASTQ files.",
    )

    parser.add_argument(
        "--reads",
        nargs="+",
        type=int,
        required=True,
        help="Maximum number of reads to sample.",
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads.",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed.",
    )

    parser.add_argument(
        "--sample",
        default=None,
        help="Process only this sample.",
    )

    parser.add_argument(
        "--sralist",
        required=True,
        help="Master list containing SRA/ENA accessions.",
    )

    parser.add_argument(
        "--parquetdir",
        required=True,
        help="Directory containing completed Parquet files.",
    )

    args = parser.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    sralist = Path(args.sralist)
    parquetdir = Path(args.parquetdir)

    indir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------
    # Read accession list
    # --------------------------------------------------------

    with open(sralist) as handle:

        accessions = [
            line.strip()
            for line in handle
            if line.strip()
        ]

    if args.sample is not None:

        if args.sample not in accessions:
            print(
                f"[WARNING] {args.sample} is not present "
                f"in {sralist}"
            )

        accessions = [args.sample]

    # --------------------------------------------------------
    # Process requested accessions
    # --------------------------------------------------------

    print("=" * 70)
    print("Checking Parquet files and raw FASTQ files")
    print("=" * 70)

    for accession in accessions:

        # ----------------------------------------------------
        # First check whether ALL requested read depths are
        # already completed.
        # ----------------------------------------------------

        all_complete = all(
            parquet_exists(
                parquetdir,
                accession,
                reads,
            )
            for reads in args.reads
        )

        if all_complete:

            print(
                f"[SKIP] {accession}: "
                f"Parquet already exists."
            )

            continue

        # ----------------------------------------------------
        # Check whether raw FASTQ already exists.
        # ----------------------------------------------------

        samples = identify_samples(indir)

        if accession in samples:

            print(
                f"[LOCAL] {accession}: "
                f"raw FASTQ already exists."
            )

        else:

            print()
            print(
                f"[DOWNLOAD] {accession}: "
                f"raw FASTQ not found."
            )

            try:

                download_sample(
                    accession,
                    indir,
                )

            except requests.exceptions.RequestException as exc:

                print(
                    f"[DOWNLOAD FAILED] {accession}"
                )
                print(f"  {exc}")

            except Exception as exc:

                print(
                    f"[DOWNLOAD FAILED] {accession}"
                )
                print(f"  {exc}")

    # --------------------------------------------------------
    # Re-scan after downloads
    # --------------------------------------------------------

    samples = identify_samples(indir)

    # --------------------------------------------------------
    # Subsample
    # --------------------------------------------------------

    for accession in accessions:

        if accession not in samples:

            # It may already have all Parquet files, in which
            # case this should not happen because of the
            # earlier skip.
            print(
                f"[MISSING RAW READS] "
                f"{accession}: no usable FASTQ found "
                f"after dnload attempt."
            )

            continue

        for reads in args.reads:

            # ------------------------------------------------
            # Parquet check
            # ------------------------------------------------

            if parquet_exists(
                parquetdir,
                accession,
                reads,
            ):

                print(
                    f"[SKIP] {accession}: "
                    f"{reads:,} reads already in Parquet."
                )

                continue

            # ------------------------------------------------
            # Check whether output already exists
            # ------------------------------------------------

            if samples[accession]["type"] == "paired":

                output_r1 = (
                    outdir
                    / f"{accession}_R1.fastq.gz"
                )

                output_r2 = (
                    outdir
                    / f"{accession}_R2.fastq.gz"
                )

                if (
                    validate_fastq(output_r1)
                    and validate_fastq(output_r2)
                ):

                    print(
                        f"[SKIP] {accession}: "
                        f"{reads:,}-read subsample "
                        f"already exists."
                    )

                    continue

            else:

                output = (
                    outdir
                    / f"{accession}.fastq.gz"
                )

                if validate_fastq(output):

                    print(
                        f"[SKIP] {accession}: "
                        f"{reads:,}-read subsample "
                        f"already exists."
                    )

                    continue

            # ------------------------------------------------
            # Perform sampling
            # ------------------------------------------------

            try:

                process_sample(
                    accession,
                    samples[accession],
                    outdir,
                    reads,
                    args.seed,
                )

            except Exception as exc:

                print(
                    f"[FAILED] {accession}: "
                    f"subsampling failed."
                )

                print(f"  {exc}")

                return 1

    print()
    print("=" * 70)
    print("Subsampler finished")
    print("=" * 70)

    return 0


if __name__ == "__main__":
    sys.exit(main())
