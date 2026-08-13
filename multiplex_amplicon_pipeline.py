#!/usr/bin/env python3
"""Process multiple paired-end amplicon libraries as one analysis set.

The workflow supports two read-joining strategies:

* ``concat``: reverse-complement R2 and concatenate it to R1 without overlap.
* ``overlap``: merge overlapping pairs with fastp and concatenate unmerged pairs.

Reads are demultiplexed using sample-specific primer/index sequences listed in a
TSV file, renamed as ``<sample_id>.<read_number>``, pooled across libraries,
filtered against a reference sequence, oriented, dereplicated, clustered, and
finally written back into per-sample FASTA files.
"""

import argparse
import csv
import logging
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s: %(message)s")


def run(cmd, check=True):
    """Run a command and log it before execution."""
    logging.info("RUN: %s", " ".join(map(str, cmd)))
    return subprocess.run(cmd, check=check)


def check_tools(mode):
    """Verify external command-line dependencies."""
    tools = ["seqkit", "vsearch", "makeblastdb", "blastn"]
    if mode == "overlap":
        tools.append("fastp")
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        sys.exit("ERROR: missing command-line tools: " + ", ".join(missing))


def safe_filename(value):
    """Convert an arbitrary sample or library label into a filesystem-safe token."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return cleaned or "unnamed"


def find_paired_libraries(r1_dir, r2_dir, r1_glob, r1_token, r2_token):
    """Find R1 files and their matching R2 files using configurable filename tokens."""
    r1_files = sorted(Path(r1_dir).glob(r1_glob))
    libraries = []
    for r1 in r1_files:
        if r1_token not in r1.name:
            logging.warning("R1 token '%s' not found in %s; skipping", r1_token, r1.name)
            continue
        r2 = Path(r2_dir) / r1.name.replace(r1_token, r2_token, 1)
        if not r2.exists():
            logging.warning("No matching R2 file for %s", r1.name)
            continue
        library_name = r1.name.split(r1_token, 1)[0]
        libraries.append({"name": library_name, "r1": str(r1), "r2": str(r2)})
    return libraries


def parse_primers_table(path):
    """Read a tab-delimited table: Library, SampleID, ..., ForwardPrimer, ReversePrimer."""
    primers = defaultdict(list)
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row or not row[0].strip():
                continue
            if len(row) < 5:
                logging.warning("Skipping primer-table row %d with fewer than 5 columns", line_number)
                continue
            library_name = row[0].strip()
            sample_id = row[1].strip()
            forward = row[3].strip()
            reverse = row[4].strip()
            if not sample_id or not forward or not reverse:
                logging.warning("Skipping incomplete primer-table row %d", line_number)
                continue
            primers[library_name].append(
                {"sample_id": sample_id, "forward": forward, "reverse": reverse}
            )
    return primers


def concatenate_files(output_file, input_files):
    """Concatenate existing files in the supplied order."""
    with open(output_file, "wb") as destination:
        for input_file in input_files:
            input_path = Path(input_file)
            if input_path.exists() and input_path.stat().st_size > 0:
                with input_path.open("rb") as source:
                    shutil.copyfileobj(source, destination)


def sample_id_from_read_name(header):
    """Recover the sample ID from a generated '<sample_id>.<number>' read name."""
    token = header.lstrip(">").strip().split()[0]
    match = re.fullmatch(r"(.+)\.(\d+)", token)
    if not match:
        raise ValueError(
            f"Read name '{token}' does not match the expected '<sample_id>.<number>' format"
        )
    return match.group(1)


def normalize_mode(mode):
    """Normalize accepted joining-mode labels."""
    aliases = {"its": "concat", "short": "overlap"}
    return aliases.get(mode, mode)


def main():
    parser = argparse.ArgumentParser(
        description="Process and jointly cluster multiple paired-end amplicon libraries.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--r1-dir", "--r1_dir", dest="r1_dir", required=True, help="Directory containing R1 FASTQ files")
    parser.add_argument("--r2-dir", "--r2_dir", dest="r2_dir", required=True, help="Directory containing R2 FASTQ files")
    parser.add_argument("--primers", required=True, help="TSV: Library, SampleID, ..., ForwardPrimer, ReversePrimer")
    parser.add_argument("--query", required=True, help="Reference sequence in FASTA format")
    parser.add_argument("--outdir", default="multi_run", help="Output directory")
    parser.add_argument(
        "--mode",
        choices=["concat", "overlap", "its", "short"],
        default="concat",
        help="Read joining mode; 'its' and 'short' are accepted aliases for compatibility",
    )
    parser.add_argument("--r1-glob", default="*_R1_*.fastq.gz", help="Glob used to find R1 files")
    parser.add_argument("--r1-token", default="_R1_", help="Token identifying R1 in paired filenames")
    parser.add_argument("--r2-token", default="_R2_", help="Token replacing --r1-token to identify R2")
    parser.add_argument("--threads", type=int, default=12, help="Number of threads")
    parser.add_argument("--fastp-min-overlap", "--fastp_min_overlap", dest="fastp_min_overlap", type=int, default=30, help="Minimum overlap for fastp merging")
    parser.add_argument("--fastp-mismatch-rate", "--fastp_mismatch_rate", dest="fastp_mismatch_rate", type=float, default=0.15, help="Maximum mismatch fraction in the overlap")
    parser.add_argument(
        "--cluster-method",
        "--cluster_method",
        dest="cluster_method",
        choices=["cluster_fast", "cluster_unoise"],
        default="cluster_fast",
        help="VSEARCH clustering method",
    )
    parser.add_argument(
        "--cluster-id",
        "--cluster_id",
        dest="cluster_id",
        type=float,
        default=1.0,
        help="Identity threshold for cluster_fast and read-to-centroid mapping",
    )
    parser.add_argument("--trim", type=int, default=3, help="Number of nucleotides trimmed from each end after target filtering")
    parser.add_argument("--blast-perc-identity", "--blast_perc_identity", dest="blast_perc_identity", type=float, default=70.0, help="Minimum BLAST percent identity")
    parser.add_argument("--blast-qcovs", "--blast_qcovs", dest="blast_qcovs", type=float, default=70.0, help="Minimum BLAST query coverage percentage")
    parser.add_argument("--blast-evalue", "--blast_evalue", dest="blast_evalue", default="1e-20", help="Maximum BLAST E-value")
    parser.add_argument("--minsize", type=int, default=5, help="Minimum abundance used by VSEARCH UNOISE")
    args = parser.parse_args()

    args.mode = normalize_mode(args.mode)

    if not 0 < args.cluster_id <= 1:
        parser.error("--cluster-id must be in the interval (0, 1]")
    if not 0 <= args.fastp_mismatch_rate <= 1:
        parser.error("--fastp-mismatch-rate must be in the interval [0, 1]")
    if args.trim < 0:
        parser.error("--trim must be non-negative")

    query = Path(args.query)
    primers_path = Path(args.primers)
    r1_dir = Path(args.r1_dir)
    r2_dir = Path(args.r2_dir)
    for path, label in [
        (query, "reference FASTA"),
        (primers_path, "primer table"),
        (r1_dir, "R1 directory"),
        (r2_dir, "R2 directory"),
    ]:
        if not path.exists():
            sys.exit(f"ERROR: {label} was not found: {path}")

    check_tools(args.mode)

    out = Path(args.outdir)
    out.mkdir(exist_ok=True, parents=True)

    logging.info("Step 1: locating paired libraries")
    libraries = find_paired_libraries(
        r1_dir, r2_dir, args.r1_glob, args.r1_token, args.r2_token
    )
    if not libraries:
        sys.exit("ERROR: no R1/R2 pairs were found.")
    logging.info("Libraries found: %s", ", ".join(lib["name"] for lib in libraries))

    logging.info("Step 2: reading the primer table")
    primer_table = parse_primers_table(primers_path)

    all_renamed_fasta = out / "all_renamed_reads.fasta"
    total_samples_processed = 0

    with all_renamed_fasta.open("w") as pooled_output:
        for index, library in enumerate(libraries, start=1):
            library_name = library["name"]
            logging.info(
                "Processing library %s (%d/%d)", library_name, index, len(libraries)
            )

            if library_name not in primer_table:
                logging.warning("No primer definitions for %s; skipping", library_name)
                continue

            temp_dir = out / f"temp_{safe_filename(library_name)}"
            temp_dir.mkdir(exist_ok=True)

            try:
                if args.mode == "concat":
                    merged_fasta = temp_dir / "concatenated_pairs.fasta"
                    r1_fasta = temp_dir / "R1.fasta"
                    r2_fasta = temp_dir / "R2.fasta"
                    r2_reverse = temp_dir / "R2_reverse_complement.fasta"
                    r1_sorted = temp_dir / "R1_sorted.fasta"
                    r2_sorted = temp_dir / "R2_sorted.fasta"

                    run(["seqkit", "fq2fa", library["r1"], "-o", str(r1_fasta)])
                    run(["seqkit", "fq2fa", library["r2"], "-o", str(r2_fasta)])
                    run(["seqkit", "seq", "-r", "-p", str(r2_fasta), "-o", str(r2_reverse)])
                    run(["seqkit", "sort", "-n", str(r1_fasta), "-o", str(r1_sorted)])
                    run(["seqkit", "sort", "-n", str(r2_reverse), "-o", str(r2_sorted)])
                    run(["seqkit", "concat", str(r1_sorted), str(r2_sorted), "-o", str(merged_fasta)])
                else:
                    logging.info("Using fastp overlap merging")
                    merged_fastq = temp_dir / "fastp_merged.fastq"
                    unmerged_r1_fastq = temp_dir / "unmerged_R1.fastq"
                    unmerged_r2_fastq = temp_dir / "unmerged_R2.fastq"
                    merged_fasta = temp_dir / "fastp_merged.fasta"

                    run(
                        [
                            "fastp",
                            "-i",
                            library["r1"],
                            "-I",
                            library["r2"],
                            "--merge",
                            "--merged_out",
                            str(merged_fastq),
                            "--out1",
                            str(unmerged_r1_fastq),
                            "--out2",
                            str(unmerged_r2_fastq),
                            "--overlap_len_require",
                            str(args.fastp_min_overlap),
                            "--overlap_diff_percent_limit",
                            str(round(args.fastp_mismatch_rate * 100)),
                            "--thread",
                            str(args.threads),
                            "--html",
                            str(temp_dir / "fastp.html"),
                            "--json",
                            str(temp_dir / "fastp.json"),
                        ]
                    )
                    run(["seqkit", "fq2fa", str(merged_fastq), "-o", str(merged_fasta)])

                    if (
                        unmerged_r1_fastq.exists()
                        and unmerged_r2_fastq.exists()
                        and unmerged_r1_fastq.stat().st_size > 0
                        and unmerged_r2_fastq.stat().st_size > 0
                    ):
                        unmerged_r1 = temp_dir / "unmerged_R1.fasta"
                        unmerged_r2 = temp_dir / "unmerged_R2.fasta"
                        unmerged_r2_reverse = temp_dir / "unmerged_R2_reverse_complement.fasta"
                        unmerged_r1_sorted = temp_dir / "unmerged_R1_sorted.fasta"
                        unmerged_r2_sorted = temp_dir / "unmerged_R2_sorted.fasta"
                        unmerged_concat = temp_dir / "unmerged_concatenated.fasta"
                        final_merged = temp_dir / "all_merged_and_unmerged.fasta"

                        run(["seqkit", "fq2fa", str(unmerged_r1_fastq), "-o", str(unmerged_r1)])
                        run(["seqkit", "fq2fa", str(unmerged_r2_fastq), "-o", str(unmerged_r2)])
                        run(["seqkit", "seq", "-r", "-p", str(unmerged_r2), "-o", str(unmerged_r2_reverse)])
                        run(["seqkit", "sort", "-n", str(unmerged_r1), "-o", str(unmerged_r1_sorted)])
                        run(["seqkit", "sort", "-n", str(unmerged_r2_reverse), "-o", str(unmerged_r2_sorted)])
                        run(["seqkit", "concat", str(unmerged_r1_sorted), str(unmerged_r2_sorted), "-o", str(unmerged_concat)])
                        concatenate_files(final_merged, [merged_fasta, unmerged_concat])
                        merged_fasta = final_merged

                for primer_set in primer_table[library_name]:
                    sample_id = primer_set["sample_id"]
                    sample_token = safe_filename(sample_id)
                    forward = primer_set["forward"]
                    reverse = primer_set["reverse"]
                    forward_hits = temp_dir / f"{sample_token}_forward.fasta"
                    sample_fasta = temp_dir / f"{sample_token}.fasta"

                    run(["seqkit", "grep", "-s", "-i", "-p", forward, str(merged_fasta), "-o", str(forward_hits)], check=False)
                    run(["seqkit", "grep", "-s", "-i", "-p", reverse, str(forward_hits), "-o", str(sample_fasta)], check=False)

                    if sample_fasta.exists() and sample_fasta.stat().st_size > 0:
                        renamed = temp_dir / f"{sample_token}_renamed.fasta"
                        run(["seqkit", "replace", "-p", "(.+)", "-r", f"{sample_id}.{{nr}}", str(sample_fasta), "-o", str(renamed)])
                        with renamed.open() as source:
                            pooled_output.write(source.read())
                        total_samples_processed += 1
                        logging.info("Added sample: %s", sample_id)
            finally:
                shutil.rmtree(temp_dir, ignore_errors=True)

    if not all_renamed_fasta.exists() or all_renamed_fasta.stat().st_size == 0:
        sys.exit(
            "ERROR: no reads passed demultiplexing. Check library names, primer definitions, and input files."
        )

    logging.info("Step 3: pooled reads from %d sample entries", total_samples_processed)

    blast_db = out / "all_reads_db"
    run(["makeblastdb", "-in", str(all_renamed_fasta), "-dbtype", "nucl", "-out", str(blast_db)])

    blast_raw = out / "blast_raw.tsv"
    run(
        [
            "blastn",
            "-query",
            str(query),
            "-db",
            str(blast_db),
            "-outfmt",
            "6 sseqid qcovs pident",
            "-max_target_seqs",
            "1000000",
            "-evalue",
            args.blast_evalue,
            "-out",
            str(blast_raw),
        ]
    )

    hit_ids = out / "blast_hits.ids"
    passed_hits = 0
    total_hits = 0
    with blast_raw.open() as source, hit_ids.open("w") as destination:
        for line in source:
            total_hits += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            seq_id, query_coverage, percent_identity = parts[0], float(parts[1]), float(parts[2])
            if percent_identity >= args.blast_perc_identity and query_coverage >= args.blast_qcovs:
                destination.write(seq_id + "\n")
                passed_hits += 1

    logging.info("BLAST hits: %d total, %d passed filters", total_hits, passed_hits)
    if passed_hits == 0:
        sys.exit("ERROR: no sequences passed the BLAST identity and query-coverage filters.")

    unique_hit_ids = out / "blast_hits.unique.ids"
    with hit_ids.open() as source:
        unique_ids = sorted({line.strip() for line in source if line.strip()})
    with unique_hit_ids.open("w") as destination:
        destination.write("\n".join(unique_ids) + "\n")

    filtered = out / "target_reads.fasta"
    run(["seqkit", "grep", "-f", str(unique_hit_ids), str(all_renamed_fasta), "-o", str(filtered)])
    if not filtered.exists() or filtered.stat().st_size == 0:
        sys.exit("ERROR: target-read extraction produced an empty FASTA file.")

    if args.trim > 0:
        trimmed = out / "trimmed.fasta"
        run(["seqkit", "subseq", "-r", f"{args.trim + 1}:-{args.trim + 1}", str(filtered), "-o", str(trimmed)])
    else:
        trimmed = filtered

    oriented = out / "oriented.fasta"
    run(["vsearch", "--orient", str(trimmed), "--db", str(query), "--fastaout", str(oriented), "--threads", str(args.threads), "--quiet"])
    if not oriented.exists() or oriented.stat().st_size == 0:
        sys.exit("ERROR: sequence orientation produced no output.")

    dereplicated = out / "dereplicated.fasta"
    run(["vsearch", "--derep_fulllength", str(oriented), "--sizeout", "--output", str(dereplicated), "--threads", str(args.threads), "--quiet"])

    if args.cluster_method == "cluster_unoise":
        centroids = out / "unoise_ASV.fasta"
        run(["vsearch", "--cluster_unoise", str(dereplicated), "--centroids", str(centroids), "--minsize", str(args.minsize), "--unoise_alpha", "4.0", "--threads", str(args.threads), "--quiet"])
    else:
        centroids = out / "sequence_variants.fasta"
        run(["vsearch", "--cluster_fast", str(dereplicated), "--id", str(args.cluster_id), "--centroids", str(centroids), "--threads", str(args.threads), "--quiet"])

    if not centroids.exists() or centroids.stat().st_size == 0:
        sys.exit("ERROR: clustering produced no centroids.")

    mapping = out / "mapping.tsv"
    run(["vsearch", "--usearch_global", str(oriented), "--db", str(centroids), "--id", str(args.cluster_id), "--blast6out", str(mapping), "--threads", str(args.threads), "--quiet"])

    matched_read_ids = out / "matched_read_ids.txt"
    with mapping.open() as source, matched_read_ids.open("w") as destination:
        for line in source:
            if line.strip():
                destination.write(line.split("\t", 1)[0] + "\n")

    final_reads = out / "final_mapped_reads.fasta"
    run(["seqkit", "grep", "-n", "-f", str(matched_read_ids), str(oriented), "-o", str(final_reads)])

    fastas_dir = out / "fastas"
    fastas_dir.mkdir(exist_ok=True)
    sample_names = set()
    with final_reads.open() as handle:
        for line in handle:
            if line.startswith(">"):
                try:
                    sample_names.add(sample_id_from_read_name(line))
                except ValueError as exc:
                    sys.exit(f"ERROR: {exc}")

    for sample in sorted(sample_names):
        sample_file = fastas_dir / f"{safe_filename(sample)}.fasta"
        pattern = rf"^{re.escape(sample)}\."
        run(["seqkit", "grep", "-r", "-p", pattern, str(final_reads), "-o", str(sample_file)], check=False)

    logging.info("Done. Results: %s", out)
    logging.info("Pooled reads: %s", all_renamed_fasta)
    logging.info("Target reads: %s", filtered)
    logging.info("Dereplicated sequences: %s", dereplicated)
    logging.info("Final centroids: %s", centroids)
    logging.info("Mapped reads: %s", final_reads)
    logging.info("Per-sample FASTA files: %s", fastas_dir)


if __name__ == "__main__":
    main()
