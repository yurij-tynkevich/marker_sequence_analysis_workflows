#!/usr/bin/env python3
"""Calculate per-sample abundance of dereplicated sequence haplotypes.

The script dereplicates a combined FASTA file with VSEARCH, retains haplotypes
meeting a global minimum copy number, maps all reads to the complete retained
haplotype set, and reports the relative abundance of selected haplotypes within
each sample.

Input read identifiers are expected to end in ``.<read_number>`` by default,
for example ``Sample_A.123``. A custom regular expression can be supplied with
``--sample-regex``; its first capture group is used as the sample identifier.
"""

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


def run_command(cmd, description):
    """Run an external command and fail with a concise diagnostic message."""
    print(f"{description} ...")
    try:
        return subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as exc:
        if exc.stderr:
            print(exc.stderr.strip(), file=sys.stderr)
        raise


def check_dependencies():
    """Verify that required command-line tools are available."""
    if shutil.which("vsearch") is None:
        sys.exit("ERROR: vsearch was not found in PATH.")


def parse_sample_id(header, sample_regex):
    """Extract a sample identifier from a FASTA header."""
    token = header.lstrip(">").strip().split()[0]
    match = re.search(sample_regex, token)
    if not match or match.lastindex is None:
        raise ValueError(
            f"Could not extract a sample ID from header '{token}' using "
            f"--sample-regex '{sample_regex}'. The expression must contain "
            "at least one capture group."
        )
    return match.group(1)


def haplotype_number(name):
    """Return the numeric part of a VSEARCH H_<n> label for natural sorting."""
    match = re.fullmatch(r"H_(\d+)", name)
    if not match:
        raise ValueError(f"Unexpected haplotype label: {name}")
    return int(match.group(1))


def parse_haplotype_selection(items):
    """Expand labels such as H_1-H_20 and return a naturally sorted list."""
    selected = set()
    for item in items:
        range_match = re.fullmatch(r"H_(\d+)-H_(\d+)", item)
        if range_match:
            start, end = map(int, range_match.groups())
            if end < start:
                raise ValueError(f"Invalid haplotype range: {item}")
            selected.update(f"H_{i}" for i in range(start, end + 1))
        else:
            haplotype_number(item)
            selected.add(item)
    return sorted(selected, key=haplotype_number)


def main():
    parser = argparse.ArgumentParser(
        description="Calculate per-sample abundance of dereplicated sequence haplotypes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("reads", help="Combined FASTA containing reads from all samples")
    parser.add_argument(
        "--sample-regex",
        default=r"^(.+)\.\d+$",
        help=(
            "Regular expression used to extract the sample ID from each read name; "
            "the first capture group is used"
        ),
    )
    parser.add_argument(
        "--min-size",
        type=int,
        default=5,
        help="Minimum global number of identical reads required for a valid haplotype",
    )
    parser.add_argument(
        "--top",
        "-t",
        type=int,
        default=20,
        help="Number of globally most abundant haplotypes to report when --haplotypes is not used",
    )
    parser.add_argument(
        "--haplotypes",
        "-hp",
        nargs="+",
        default=[],
        help="Haplotypes to report explicitly, e.g. H_1-H_20 H_25",
    )
    parser.add_argument(
        "--identity",
        "-i",
        type=float,
        default=0.97,
        help="VSEARCH identity threshold for mapping reads to valid haplotypes",
    )
    parser.add_argument("--threads", type=int, default=8, help="Number of VSEARCH threads")
    parser.add_argument(
        "--output",
        "-o",
        default=None,
        help="Output TSV/XLSX table; generated automatically when omitted",
    )
    parser.add_argument(
        "--haps-out",
        default=None,
        help="Output FASTA containing all valid haplotypes; generated automatically when omitted",
    )
    args = parser.parse_args()

    if args.min_size < 1:
        parser.error("--min-size must be at least 1")
    if not 0 < args.identity <= 1:
        parser.error("--identity must be in the interval (0, 1]")
    if args.top < 1:
        parser.error("--top must be at least 1")

    reads = Path(args.reads)
    if not reads.is_file():
        sys.exit(f"ERROR: input FASTA was not found: {reads}")

    check_dependencies()

    output_path = Path(args.output or f"haplotype_table_min{args.min_size}.tsv")
    haplotypes_path = Path(args.haps_out or f"haplotypes_min{args.min_size}.fasta")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    haplotypes_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="haplotype_distribution_") as tmpdir:
        tmpdir = Path(tmpdir)
        derep_fasta = tmpdir / "dereplicated.fasta"
        sample_reads = tmpdir / "sample_labeled_reads.fasta"
        uc_file = tmpdir / "mapping.uc"

        run_command(
            [
                "vsearch",
                "--derep_fulllength",
                str(reads),
                "--sizeout",
                "--relabel",
                "H_",
                "--output",
                str(derep_fasta),
                "--threads",
                str(args.threads),
            ],
            "Dereplicating reads",
        )

        valid_count = 0
        with derep_fasta.open() as source, haplotypes_path.open("w") as destination:
            keep = False
            for line in source:
                if line.startswith(">"):
                    match = re.search(r";size=(\d+);?", line)
                    if not match:
                        raise RuntimeError(f"Could not parse VSEARCH size annotation: {line.strip()}")
                    keep = int(match.group(1)) >= args.min_size
                    if keep:
                        valid_count += 1
                if keep:
                    destination.write(line)

        if valid_count == 0:
            sys.exit(
                f"ERROR: no haplotypes met the global minimum size of {args.min_size}."
            )

        read_counter = 0
        with reads.open() as source, sample_reads.open("w") as destination:
            for line in source:
                if line.startswith(">"):
                    try:
                        sample = parse_sample_id(line, args.sample_regex)
                    except ValueError as exc:
                        sys.exit(f"ERROR: {exc}")
                    read_counter += 1
                    destination.write(f">{sample}|read_{read_counter}\n")
                else:
                    destination.write(line)

        run_command(
            [
                "vsearch",
                "--usearch_global",
                str(sample_reads),
                "--db",
                str(haplotypes_path),
                "--id",
                str(args.identity),
                "--uc",
                str(uc_file),
                "--threads",
                str(args.threads),
                "--maxaccepts",
                "1",
                "--maxrejects",
                "32",
                "--strand",
                "plus",
                "--quiet",
            ],
            "Mapping reads to the complete valid haplotype set",
        )

        valid_haplotypes = []
        with haplotypes_path.open() as handle:
            for line in handle:
                if line.startswith(">"):
                    valid_haplotypes.append(line[1:].split(";")[0].strip())
        valid_set = set(valid_haplotypes)

        sample_haplotype_counts = defaultdict(Counter)
        sample_mapped_totals = Counter()

        with uc_file.open() as handle:
            for line in handle:
                if not line.startswith("H"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    continue
                query_label = fields[8]
                target = fields[9].split(";")[0]
                sample = query_label.split("|read_", 1)[0]
                if target in valid_set:
                    sample_mapped_totals[sample] += 1
                    sample_haplotype_counts[sample][target] += 1

        if not sample_mapped_totals:
            sys.exit("ERROR: no reads mapped to the valid haplotype set.")

        if args.haplotypes:
            try:
                selected = parse_haplotype_selection(args.haplotypes)
            except ValueError as exc:
                parser.error(str(exc))
            missing = [hap for hap in selected if hap not in valid_set]
            if missing:
                print(
                    "WARNING: requested haplotypes not present in the valid set: "
                    + ", ".join(missing),
                    file=sys.stderr,
                )
        else:
            global_counts = Counter()
            for counts in sample_haplotype_counts.values():
                global_counts.update(counts)
            selected = [hap for hap, _ in global_counts.most_common(args.top)]
            selected.sort(key=haplotype_number)

        rows = []
        for sample in sorted(sample_mapped_totals):
            total = sample_mapped_totals[sample]
            row = {"Sample": sample, "Mapped_reads": total}
            selected_total = 0
            for haplotype in selected:
                count = sample_haplotype_counts[sample].get(haplotype, 0)
                row[f"hap_{haplotype}"] = round(count * 100 / total, 1)
                selected_total += count
            row["Other"] = round(max(total - selected_total, 0) * 100 / total, 1)
            rows.append(row)

        columns = ["Sample", "Mapped_reads"] + [f"hap_{hap}" for hap in selected] + ["Other"]
        dataframe = pd.DataFrame(rows).reindex(columns=columns, fill_value=0.0)

        if output_path.suffix.lower() in {".xlsx", ".xls"}:
            dataframe.to_excel(output_path, index=False)
        else:
            dataframe.to_csv(output_path, sep="\t", index=False)

    print("Done.")
    print(f"Valid haplotypes (global size >= {args.min_size}): {haplotypes_path}")
    print(f"Mapped abundance table: {output_path}")
    print(dataframe.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
