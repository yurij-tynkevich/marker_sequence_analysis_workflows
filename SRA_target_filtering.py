#!/usr/bin/env python3
"""Extract target-matching read segments from local FASTQ files or an SRA run.

Short-read libraries are mapped with BWA-MEM2. Long-read libraries are mapped
with minimap2 using platform-specific presets for Oxford Nanopore, PacBio HiFi,
or PacBio CLR data. Only the aligned segment of each mapped read is written to
the merged FASTQ output.

For paired-end data, R1 and R2 are retained independently. Deduplication is
performed within each input stream rather than globally by query name, so a
mapped R2 mate is not discarded merely because its R1 mate has the same name.
"""

import argparse
import logging
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s: %(message)s")


def run(cmd, check=True, shell=False):
    """Run an external command and return the completed process."""
    if shell:
        logging.info("RUN: %s", cmd)
    else:
        logging.info("RUN: %s", " ".join(map(str, cmd)))
    return subprocess.run(cmd, check=check, shell=shell, executable="/bin/bash" if shell else None)


def check_tools(mode, srr_present, use_fastp):
    """Verify required command-line tools."""
    tools = ["samtools", "seqtk", "seqkit"]
    if use_fastp:
        tools.append("fastp")
    if srr_present:
        tools.extend(["prefetch", "fastq-dump"])
    tools.append("bwa-mem2" if mode == "short" else "minimap2")
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        sys.exit("ERROR: missing command-line tools: " + ", ".join(missing))


def download_sra(srr, out_dir, paired=False):
    """Download an SRA run and convert it to gzipped FASTQ."""
    download_dir = out_dir / "raw_sra"
    download_dir.mkdir(exist_ok=True, parents=True)

    logging.info("Downloading %s from SRA", srr)
    run(["prefetch", srr])

    command = ["fastq-dump", "--gzip", srr, "--outdir", str(download_dir)]
    if paired:
        command.append("--split-files")
    run(command)

    raw_files = sorted(download_dir.glob(f"{srr}*.fastq.gz"))
    if not raw_files:
        sys.exit(f"ERROR: no FASTQ files were produced for {srr}.")
    if paired and len(raw_files) != 2:
        sys.exit(
            f"ERROR: --paired was requested, but {len(raw_files)} FASTQ files were produced for {srr}."
        )
    return raw_files


def sample_fastq_files(input_files, out_dir, max_reads, paired=False, seed=100):
    """Randomly subsample FASTQ files with seqtk using a fixed seed."""
    if not max_reads:
        return list(input_files)

    logging.info("Randomly sampling %d reads per input stream", max_reads)
    sampled_files = []
    for index, input_file in enumerate(input_files, start=1):
        suffix = f"R{index}" if paired else f"stream{index}"
        sampled = out_dir / f"sampled_{suffix}.fastq.gz"
        command = (
            f"seqtk sample -s{seed} {shlex.quote(str(input_file))} {max_reads} "
            f"| gzip > {shlex.quote(str(sampled))}"
        )
        run(command, shell=True)
        sampled_files.append(sampled)
    return sampled_files


def preprocess_fastq(input_files, out_dir, paired=False, threads=8):
    """Quality-filter FASTQ input with fastp."""
    out_dir.mkdir(exist_ok=True, parents=True)

    if paired:
        if len(input_files) != 2:
            raise ValueError("Paired-end preprocessing requires exactly two FASTQ files.")
        in1, in2 = input_files
        out1 = out_dir / "clean_R1.fastq.gz"
        out2 = out_dir / "clean_R2.fastq.gz"
        run(
            [
                "fastp",
                "-i",
                str(in1),
                "-I",
                str(in2),
                "-o",
                str(out1),
                "-O",
                str(out2),
                "-q",
                "20",
                "-l",
                "50",
                "--thread",
                str(threads),
                "--html",
                str(out_dir / "fastp.html"),
                "--json",
                str(out_dir / "fastp.json"),
            ]
        )
        return [out1, out2]

    clean_files = []
    for index, input_file in enumerate(input_files, start=1):
        out_file = out_dir / f"clean_stream{index}.fastq.gz"
        run(
            [
                "fastp",
                "-i",
                str(input_file),
                "-o",
                str(out_file),
                "-q",
                "20",
                "-l",
                "50",
                "--thread",
                str(threads),
                "--html",
                str(out_dir / f"fastp_stream{index}.html"),
                "--json",
                str(out_dir / f"fastp_stream{index}.json"),
            ]
        )
        clean_files.append(out_file)
    return clean_files


def extract_aligned_segments(
    bam_path,
    output_fastq,
    source_id,
    record_suffix="",
    min_aligned_length=50,
    min_mapq=0,
):
    """Append one aligned segment per query from a BAM file to a merged FASTQ."""
    try:
        import pysam
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
    except ImportError as exc:
        sys.exit(
            f"ERROR: missing Python dependency: {exc.name}. "
            "Install requirements.txt before running this workflow."
        )

    seen_queries = set()
    extracted = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam, output_fastq.open("a") as handle:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.query_sequence is None:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if read.query_name in seen_queries:
                continue

            start = read.query_alignment_start
            end = read.query_alignment_end
            if start is None or end is None or (end - start) < min_aligned_length:
                continue

            sequence = read.query_sequence[start:end]
            qualities = (
                read.query_qualities[start:end]
                if read.query_qualities is not None
                else [30] * len(sequence)
            )
            base_query_name = re.sub(r"/[12]$", "", read.query_name) if record_suffix else read.query_name
            record_id = f"{base_query_name}{record_suffix}"
            description = f"source={source_id} MAPQ={read.mapping_quality}"
            record = SeqRecord(
                Seq(sequence),
                id=record_id,
                description=description,
                letter_annotations={"phred_quality": qualities},
            )
            SeqIO.write(record, handle, "fastq")
            seen_queries.add(read.query_name)
            extracted += 1

    return extracted


def minimap_preset(long_read_type):
    """Translate a platform label into the corresponding minimap2 preset."""
    return {
        "ont": "map-ont",
        "pacbio-hifi": "map-hifi",
        "pacbio-clr": "map-pb",
    }[long_read_type]


def main():
    parser = argparse.ArgumentParser(
        description="Extract target-matching read segments from local FASTQ files or SRA.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--haplotypes", required=True, help="Target sequences in FASTA format")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--srr", help="SRA run accession")
    source.add_argument("--input", nargs="+", help="Local FASTQ or FASTQ.GZ input files")
    parser.add_argument("--paired", action="store_true", help="Treat input as paired-end short reads")
    parser.add_argument("--skip-fastp", action="store_true", help="Skip fastp preprocessing")
    parser.add_argument("--max-reads", type=int, help="Randomly sample this many reads from each input stream")
    parser.add_argument("--sample-seed", type=int, default=100, help="Random seed passed to seqtk sample")
    parser.add_argument("--mode", choices=["short", "long"], default="short", help="Mapper family")
    parser.add_argument(
        "--long-read-type",
        choices=["ont", "pacbio-hifi", "pacbio-clr"],
        default="ont",
        help="Long-read platform used to select the minimap2 preset",
    )
    parser.add_argument("--threads", type=int, default=12, help="Number of mapper/fastp threads")
    parser.add_argument("--bwa-score-threshold", type=int, default=30, help="BWA-MEM2 minimum alignment score (-T)")
    parser.add_argument("--min-aligned-length", type=int, default=50, help="Minimum aligned query length retained in the output")
    parser.add_argument("--min-mapq", type=int, default=0, help="Minimum mapping quality retained in the output")
    parser.add_argument("--outdir", default="results_target", help="Output directory")
    args = parser.parse_args()

    if args.max_reads is not None and args.max_reads < 1:
        parser.error("--max-reads must be at least 1")
    if args.min_aligned_length < 1:
        parser.error("--min-aligned-length must be at least 1")
    if args.mode == "long" and args.paired:
        parser.error("--paired is supported only with --mode short")

    target_fasta = Path(args.haplotypes)
    if not target_fasta.is_file():
        sys.exit(f"ERROR: target FASTA was not found: {target_fasta}")

    local_inputs = [Path(path) for path in args.input] if args.input else None
    if local_inputs:
        missing = [str(path) for path in local_inputs if not path.is_file()]
        if missing:
            sys.exit("ERROR: input FASTQ files were not found: " + ", ".join(missing))
        if args.paired and len(local_inputs) != 2:
            parser.error("--paired requires exactly two local FASTQ files")

    use_fastp = not args.skip_fastp
    check_tools(args.mode, bool(args.srr), use_fastp)

    out = Path(args.outdir)
    out.mkdir(exist_ok=True, parents=True)
    merged_fastq = out / "merged_targets.fastq"
    merged_fastq.write_text("")

    downloaded_files = []
    if args.srr:
        current_files = download_sra(args.srr, out, paired=args.paired)
        downloaded_files = list(current_files)
    else:
        current_files = local_inputs

    sampled_files = sample_fastq_files(
        current_files,
        out,
        args.max_reads,
        paired=args.paired,
        seed=args.sample_seed,
    )

    if use_fastp:
        clean_dir = out / "clean"
        clean_files = preprocess_fastq(
            sampled_files,
            clean_dir,
            paired=args.paired,
            threads=args.threads,
        )
    else:
        clean_dir = None
        clean_files = list(sampled_files)

    filter_tmp = out / "filter_tmp"
    shutil.rmtree(filter_tmp, ignore_errors=True)
    filter_tmp.mkdir(exist_ok=True)

    reference_copy = filter_tmp / "targets.fasta"
    shutil.copyfile(target_fasta, reference_copy)

    if args.mode == "short":
        run(["bwa-mem2", "index", str(reference_copy)])
        mapping_reference = reference_copy
    else:
        mapping_reference = filter_tmp / "targets.mmi"
        run(["minimap2", "-d", str(mapping_reference), str(reference_copy)])

    total_found = 0
    for index, fastq_file in enumerate(clean_files, start=1):
        source_id = f"R{index}" if args.paired else f"stream{index}"
        suffix = f"/{index}" if args.paired else ""
        bam = filter_tmp / f"mapped_{index}.bam"
        logging.info("Mapping %s", fastq_file.name)

        if args.mode == "short":
            command = (
                f"bwa-mem2 mem -t {args.threads} -a -T {args.bwa_score_threshold} "
                f"{shlex.quote(str(mapping_reference))} {shlex.quote(str(fastq_file))} "
                f"| samtools view -b - "
                f"| samtools sort -o {shlex.quote(str(bam))} -"
            )
        else:
            preset = minimap_preset(args.long_read_type)
            command = (
                f"minimap2 -ax {preset} -t {args.threads} "
                f"{shlex.quote(str(mapping_reference))} {shlex.quote(str(fastq_file))} "
                f"| samtools view -b - "
                f"| samtools sort -o {shlex.quote(str(bam))} -"
            )
        run(command, shell=True)

        count = extract_aligned_segments(
            bam,
            merged_fastq,
            source_id=source_id,
            record_suffix=suffix,
            min_aligned_length=args.min_aligned_length,
            min_mapq=args.min_mapq,
        )
        total_found += count
        logging.info("Retained %d target reads from %s", count, source_id)

    shutil.rmtree(filter_tmp, ignore_errors=True)

    if args.srr and sampled_files != downloaded_files:
        for path in downloaded_files:
            path.unlink(missing_ok=True)
        raw_dir = out / "raw_sra"
        if raw_dir.exists() and not any(raw_dir.iterdir()):
            raw_dir.rmdir()

    logging.info("Done")
    logging.info("Unique target reads retained across input streams: %d", total_found)
    if clean_dir is not None:
        logging.info("Preprocessed reads: %s", clean_dir)
    logging.info("Merged target reads: %s", merged_fastq)
    run(["seqkit", "stats", str(merged_fastq)])


if __name__ == "__main__":
    main()
