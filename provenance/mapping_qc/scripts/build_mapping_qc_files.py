#!/usr/bin/env python3

import argparse
import csv
import gzip
import hashlib
import json
import re
import shlex
import statistics
import subprocess
import sys
from pathlib import Path


def run_command(command):
    return subprocess.check_output(
        command,
        text=True,
        stderr=subprocess.PIPE,
    ).strip()


def normalize_sample_name(value):
    name = Path(value).name

    for suffix in (
        ".sorted.bam",
        ".bam",
        ".sam",
        ".fastq.gz",
        ".fq.gz",
        ".fastq",
        ".fq",
    ):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break

    return name


def format_percent(numerator, denominator):
    if denominator == 0:
        return "NA"
    return f"{100.0 * numerator / denominator:.4f}"


def write_tsv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def count_fastq_reads(path):
    if path.name.endswith(".gz"):
        command = (
            f"gzip -cd -- {shlex.quote(str(path))} "
            "| wc -l"
        )
        lines = int(
            subprocess.check_output(
                ["bash", "-lc", command],
                text=True,
            ).strip()
        )
    else:
        lines = int(run_command(["wc", "-l", str(path)]).split()[0])

    if lines % 4 != 0:
        raise RuntimeError(
            f"FASTQ line count is not divisible by four: {path} ({lines})"
        )

    return lines // 4


def find_fastp_json(json_files, sample):
    search_names = [
        sample,
        sample.removesuffix("_RAW"),
    ]

    candidates = [
        path
        for path in json_files
        if any(name in path.name for name in search_names)
    ]

    valid = []

    for path in candidates:
        try:
            with path.open(encoding="utf-8") as handle:
                data = json.load(handle)

            before = data["summary"]["before_filtering"]["total_reads"]
            after = data["summary"]["after_filtering"]["total_reads"]

            valid.append((path, int(before), int(after)))
        except (KeyError, json.JSONDecodeError, OSError, TypeError):
            continue

    if not valid:
        raise FileNotFoundError(
            f"No valid fastp JSON found for sample {sample}"
        )

    valid.sort(
        key=lambda item: (
            0 if sample in item[0].name else 1,
            len(str(item[0])),
        )
    )

    return valid[0]


def find_host_fastq(run_dir, sample):
    direct_candidates = [
        run_dir / "temp" / f"{sample}_host.fq",
        run_dir / "temp" / f"{sample}_host.fastq",
        run_dir / "temp" / f"{sample}_host.fq.gz",
        run_dir / "temp" / f"{sample}_host.fastq.gz",
    ]

    for path in direct_candidates:
        if path.is_file():
            return path

    matches = []

    for suffix in (
        f"{sample}_host.fq",
        f"{sample}_host.fastq",
        f"{sample}_host.fq.gz",
        f"{sample}_host.fastq.gz",
    ):
        matches.extend(run_dir.rglob(suffix))

    matches = sorted(set(matches), key=lambda path: len(str(path)))

    if not matches:
        raise FileNotFoundError(
            f"No host-associated FASTQ found for sample {sample}"
        )

    return matches[0]


def extract_mapping_metrics(run_dir, output_dir):
    bam_files = sorted(
        (run_dir / "temp" / "BAM").glob("*.sorted.bam")
    )

    if not bam_files:
        raise FileNotFoundError(
            f"No sorted BAM files found under {run_dir / 'temp' / 'BAM'}"
        )

    json_files = list(run_dir.rglob("*.json"))
    rows = []

    for index, bam in enumerate(bam_files, start=1):
        sample = normalize_sample_name(bam.name)

        print(
            f"[Mapping QC] [{index}/{len(bam_files)}] {sample}",
            file=sys.stderr,
        )

        json_path, raw_reads, postfastp_reads = find_fastp_json(
            json_files,
            sample,
        )

        host_fastq = find_host_fastq(run_dir, sample)
        host_reads = count_fastq_reads(host_fastq)

        primary_mapped = int(
            run_command(
                [
                    "samtools",
                    "view",
                    "-c",
                    "-F",
                    "2308",
                    str(bam),
                ]
            )
        )

        tissue = (
            "Liver"
            if "_LIVER_" in sample
            else "Telencephalon"
            if "_TEL_" in sample
            else "Unknown"
        )

        rows.append(
            {
                "sample": sample,
                "tissue": tissue,
                "raw_reads_before_fastp": raw_reads,
                "reads_after_fastp": postfastp_reads,
                "host_classified_reads": host_reads,
                "primary_mapped_reads": primary_mapped,
                "host_partition_percent": (
                    100.0 * host_reads / postfastp_reads
                    if postfastp_reads
                    else 0.0
                ),
                "mapping_percent_of_host_reads": (
                    100.0 * primary_mapped / host_reads
                    if host_reads
                    else 0.0
                ),
                "mapping_percent_of_postfastp_reads": (
                    100.0 * primary_mapped / postfastp_reads
                    if postfastp_reads
                    else 0.0
                ),
                "fastp_json": str(json_path),
                "host_fastq": str(host_fastq),
                "bam": str(bam),
            }
        )

    table_rows = []

    for row in rows:
        table_rows.append(
            [
                row["sample"],
                row["tissue"],
                row["raw_reads_before_fastp"],
                row["reads_after_fastp"],
                row["host_classified_reads"],
                row["primary_mapped_reads"],
                f'{row["host_partition_percent"]:.4f}',
                f'{row["mapping_percent_of_host_reads"]:.4f}',
                f'{row["mapping_percent_of_postfastp_reads"]:.4f}',
                "PASS",
            ]
        )

    write_tsv(
        output_dir / "host_mapping_efficiency.tsv",
        [
            "sample",
            "tissue",
            "raw_reads_before_fastp",
            "reads_after_fastp",
            "host_classified_reads",
            "primary_mapped_reads",
            "host_partition_percent",
            "mapping_percent_of_host_reads",
            "mapping_percent_of_postfastp_reads",
            "bam_status",
        ],
        table_rows,
    )

    total_raw = sum(row["raw_reads_before_fastp"] for row in rows)
    total_post = sum(row["reads_after_fastp"] for row in rows)
    total_host = sum(row["host_classified_reads"] for row in rows)
    total_mapped = sum(row["primary_mapped_reads"] for row in rows)

    summary_rows = [
        ["libraries", len(rows)],
        ["raw_reads_before_fastp", total_raw],
        ["reads_after_fastp", total_post],
        ["host_classified_reads", total_host],
        ["primary_mapped_reads", total_mapped],
        [
            "pooled_host_partition_percent",
            format_percent(total_host, total_post),
        ],
        [
            "pooled_mapping_percent_of_host_reads",
            format_percent(total_mapped, total_host),
        ],
        [
            "pooled_mapping_percent_of_postfastp_reads",
            format_percent(total_mapped, total_post),
        ],
    ]

    for metric in (
        "host_partition_percent",
        "mapping_percent_of_host_reads",
        "mapping_percent_of_postfastp_reads",
    ):
        values = [row[metric] for row in rows]

        summary_rows.extend(
            [
                [f"{metric}_mean", f"{statistics.mean(values):.4f}"],
                [f"{metric}_median", f"{statistics.median(values):.4f}"],
                [f"{metric}_minimum", f"{min(values):.4f}"],
                [f"{metric}_maximum", f"{max(values):.4f}"],
            ]
        )

    write_tsv(
        output_dir / "host_mapping_summary.tsv",
        ["metric", "value"],
        summary_rows,
    )

    tissue_rows = []

    for tissue in ("Liver", "Telencephalon"):
        subset = [row for row in rows if row["tissue"] == tissue]

        if not subset:
            continue

        for metric in (
            "host_partition_percent",
            "mapping_percent_of_host_reads",
            "mapping_percent_of_postfastp_reads",
        ):
            values = [row[metric] for row in subset]

            tissue_rows.append(
                [
                    tissue,
                    metric,
                    len(values),
                    f"{statistics.mean(values):.4f}",
                    f"{statistics.median(values):.4f}",
                    f"{min(values):.4f}",
                    f"{max(values):.4f}",
                ]
            )

    write_tsv(
        output_dir / "host_mapping_summary_by_tissue.tsv",
        [
            "tissue",
            "metric",
            "libraries",
            "mean",
            "median",
            "minimum",
            "maximum",
        ],
        tissue_rows,
    )

    return rows


def parse_featurecounts_summary(run_dir, mapping_rows, output_dir):
    matches = sorted(run_dir.rglob("host_counts.txt.summary"))

    if not matches:
        raise FileNotFoundError(
            "Could not find host_counts.txt.summary"
        )

    summary_path = matches[0]

    with summary_path.open(encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        sample_columns = header[1:]
        status_rows = {}

        for row in reader:
            if not row:
                continue

            status_rows[row[0]] = [
                int(float(value))
                for value in row[1:]
            ]

    assigned_values = status_rows.get("Assigned")

    if assigned_values is None:
        raise RuntimeError(
            f"No Assigned row found in {summary_path}"
        )

    mapping_lookup = {
        row["sample"]: row
        for row in mapping_rows
    }

    output_rows = []
    assigned_lookup = {}

    for index, column in enumerate(sample_columns):
        sample = normalize_sample_name(column)
        total_records = sum(
            values[index]
            for values in status_rows.values()
        )
        assigned = assigned_values[index]
        primary_mapped = mapping_lookup.get(
            sample,
            {},
        ).get("primary_mapped_reads", 0)

        assigned_lookup[sample] = assigned

        output_rows.append(
            [
                sample,
                total_records,
                assigned,
                format_percent(assigned, total_records),
                primary_mapped,
                format_percent(assigned, primary_mapped),
            ]
        )

    write_tsv(
        output_dir / "featurecounts_assignment_summary.tsv",
        [
            "sample",
            "featurecounts_input_records",
            "featurecounts_assigned_reads",
            "assigned_percent_of_featurecounts_records",
            "primary_mapped_reads",
            "assigned_percent_of_primary_mapped_reads",
        ],
        output_rows,
    )

    return assigned_lookup, summary_path


def parse_gtf_gene_ids(gtf_path):
    gene_pattern = re.compile(r'gene_id\s+"([^"]+)"')
    gene_ids = set()

    opener = gzip.open if gtf_path.name.endswith(".gz") else open

    with opener(gtf_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                continue

            match = gene_pattern.search(fields[8])

            if match:
                gene_ids.add(match.group(1))

    return gene_ids


def parse_count_matrix(
    run_dir,
    gtf_path,
    mapping_rows,
    assigned_lookup,
    output_dir,
):
    matches = sorted(
        path
        for path in run_dir.rglob("host_counts.txt")
        if path.is_file()
    )

    if not matches:
        raise FileNotFoundError("Could not find host_counts.txt")

    counts_path = matches[0]

    with counts_path.open(encoding="utf-8") as handle:
        content = (
            line
            for line in handle
            if not line.startswith("#")
        )
        reader = csv.reader(content, delimiter="\t")
        header = next(reader)

        if len(header) < 7:
            raise RuntimeError(
                f"Unexpected featureCounts matrix format: {counts_path}"
            )

        sample_columns = header[6:]
        samples = [
            normalize_sample_name(column)
            for column in sample_columns
        ]

        nonzero_by_sample = {
            sample: 0
            for sample in samples
        }

        genes_in_matrix = 0
        genes_detected_any = 0

        for row in reader:
            if len(row) < 6 + len(samples):
                continue

            genes_in_matrix += 1
            counts = [
                float(value)
                for value in row[6 : 6 + len(samples)]
            ]

            detected_any = False

            for sample, value in zip(samples, counts):
                if value > 0:
                    nonzero_by_sample[sample] += 1
                    detected_any = True

            if detected_any:
                genes_detected_any += 1

    annotated_gene_ids = parse_gtf_gene_ids(gtf_path)
    annotated_genes = len(annotated_gene_ids)

    mapping_lookup = {
        row["sample"]: row
        for row in mapping_rows
    }

    rows = []

    for sample in samples:
        detected = nonzero_by_sample[sample]
        primary_mapped = mapping_lookup.get(
            sample,
            {},
        ).get("primary_mapped_reads", 0)
        assigned = assigned_lookup.get(sample, 0)

        rows.append(
            [
                sample,
                annotated_genes,
                genes_in_matrix,
                detected,
                format_percent(detected, annotated_genes),
                primary_mapped,
                assigned,
                format_percent(assigned, primary_mapped),
            ]
        )

    write_tsv(
        output_dir / "annotation_coverage_by_library.tsv",
        [
            "sample",
            "annotated_genes_in_gtf",
            "genes_in_featurecounts_matrix",
            "genes_with_nonzero_counts",
            "detected_gene_percent_of_annotated_genes",
            "primary_mapped_reads",
            "featurecounts_assigned_reads",
            "assigned_percent_of_primary_mapped_reads",
        ],
        rows,
    )

    overall_rows = [
        ["annotated_genes_in_gtf", annotated_genes],
        ["genes_in_featurecounts_matrix", genes_in_matrix],
        [
            "genes_detected_in_at_least_one_library",
            genes_detected_any,
        ],
        [
            "detected_gene_percent_of_annotated_genes",
            format_percent(genes_detected_any, annotated_genes),
        ],
        ["gtf_file", str(gtf_path)],
        ["featurecounts_matrix", str(counts_path)],
    ]

    write_tsv(
        output_dir / "annotation_coverage_summary.tsv",
        ["metric", "value"],
        overall_rows,
    )


def clean_mismatch_table(run_dir, output_dir):
    source = (
        run_dir
        / "provenance"
        / "mapping_qc"
        / "host_alignment_error_metrics.tsv"
    )

    if not source.is_file():
        print(
            f"[WARNING] Mismatch source not found: {source}",
            file=sys.stderr,
        )
        return

    with source.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    cleaned_rows = []

    for row in rows:
        cleaned_rows.append(
            [
                row["sample"],
                row["primary_mapped_reads"],
                row["reads_with_NM"],
                row["nm_tag_coverage_percent"],
                row["aligned_comparison_bases"],
                row["estimated_substitution_mismatches"],
                row["substitution_mismatch_percent"],
            ]
        )

    write_tsv(
        output_dir / "host_substitution_mismatch_by_library.tsv",
        [
            "sample",
            "primary_mapped_reads",
            "reads_with_NM",
            "nm_tag_coverage_percent",
            "aligned_comparison_bases",
            "estimated_substitution_mismatches",
            "substitution_mismatch_percent",
        ],
        cleaned_rows,
    )

    mismatch_values = [
        float(row["substitution_mismatch_percent"])
        for row in rows
    ]

    total_comparison = sum(
        int(row["aligned_comparison_bases"])
        for row in rows
    )

    total_substitutions = sum(
        int(row["estimated_substitution_mismatches"])
        for row in rows
    )

    summary_rows = [
        ["libraries", len(rows)],
        [
            "nm_tag_coverage_percent_minimum",
            f'{min(float(row["nm_tag_coverage_percent"]) for row in rows):.4f}',
        ],
        [
            "nm_tag_coverage_percent_maximum",
            f'{max(float(row["nm_tag_coverage_percent"]) for row in rows):.4f}',
        ],
        [
            "substitution_mismatch_percent_mean",
            f"{statistics.mean(mismatch_values):.4f}",
        ],
        [
            "substitution_mismatch_percent_median",
            f"{statistics.median(mismatch_values):.4f}",
        ],
        [
            "substitution_mismatch_percent_minimum",
            f"{min(mismatch_values):.4f}",
        ],
        [
            "substitution_mismatch_percent_maximum",
            f"{max(mismatch_values):.4f}",
        ],
        [
            "pooled_base_weighted_substitution_mismatch_percent",
            format_percent(total_substitutions, total_comparison),
        ],
        [
            "calculation",
            "NM minus CIGAR insertion and deletion lengths",
        ],
        [
            "scope",
            "Primary Magic-BLAST alignments against Myotis lucifugus Myoluc2.0",
        ],
    ]

    write_tsv(
        output_dir / "host_substitution_mismatch_summary.tsv",
        ["metric", "value"],
        summary_rows,
    )


def read_bam_dictionary(bam_path):
    header = run_command(
        ["samtools", "view", "-H", str(bam_path)]
    )

    sequences = {}

    for line in header.splitlines():
        if not line.startswith("@SQ\t"):
            continue

        sequence_name = None
        sequence_length = None

        for field in line.split("\t")[1:]:
            if field.startswith("SN:"):
                sequence_name = field[3:]
            elif field.startswith("LN:"):
                sequence_length = int(field[3:])

        if sequence_name is not None and sequence_length is not None:
            sequences[sequence_name] = sequence_length

    return sequences


def read_fasta_dictionary(fasta_path):
    opener = gzip.open if fasta_path.name.endswith(".gz") else open
    sequences = {}

    current_name = None
    current_length = 0

    with opener(
        fasta_path,
        "rt",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        for line in handle:
            if line.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = current_length

                current_name = line[1:].strip().split()[0]
                current_length = 0
            else:
                current_length += len(line.strip())

    if current_name is not None:
        sequences[current_name] = current_length

    return sequences


def validate_reference_dictionary(
    run_dir,
    fasta_path,
    reference_output_dir,
):
    bam_files = sorted(
        (run_dir / "temp" / "BAM").glob("*.sorted.bam")
    )

    first_bam = bam_files[0]

    print(
        "[Reference QC] Reading BAM sequence dictionary",
        file=sys.stderr,
    )
    bam_dictionary = read_bam_dictionary(first_bam)

    print(
        "[Reference QC] Reading FASTA sequence dictionary",
        file=sys.stderr,
    )
    fasta_dictionary = read_fasta_dictionary(fasta_path)

    bam_ids = set(bam_dictionary)
    fasta_ids = set(fasta_dictionary)

    only_bam = sorted(bam_ids - fasta_ids)
    only_fasta = sorted(fasta_ids - bam_ids)

    length_mismatches = sorted(
        sequence_id
        for sequence_id in bam_ids & fasta_ids
        if bam_dictionary[sequence_id]
        != fasta_dictionary[sequence_id]
    )

    exact_match = (
        not only_bam
        and not only_fasta
        and not length_mismatches
    )

    summary_rows = [
        ["bam_used", str(first_bam)],
        ["fasta_used", str(fasta_path)],
        ["bam_reference_sequences", len(bam_dictionary)],
        ["fasta_reference_sequences", len(fasta_dictionary)],
        ["sequence_ids_only_in_bam", len(only_bam)],
        ["sequence_ids_only_in_fasta", len(only_fasta)],
        ["sequence_length_mismatches", len(length_mismatches)],
        ["dictionary_exact_match", "PASS" if exact_match else "FAIL"],
    ]

    write_tsv(
        reference_output_dir / "reference_dictionary_validation.tsv",
        ["metric", "value"],
        summary_rows,
    )

    mismatch_rows = []

    for sequence_id in only_bam:
        mismatch_rows.append(
            [
                "only_in_bam",
                sequence_id,
                bam_dictionary[sequence_id],
                "NA",
            ]
        )

    for sequence_id in only_fasta:
        mismatch_rows.append(
            [
                "only_in_fasta",
                sequence_id,
                "NA",
                fasta_dictionary[sequence_id],
            ]
        )

    for sequence_id in length_mismatches:
        mismatch_rows.append(
            [
                "length_mismatch",
                sequence_id,
                bam_dictionary[sequence_id],
                fasta_dictionary[sequence_id],
            ]
        )

    write_tsv(
        reference_output_dir / "reference_dictionary_mismatches.tsv",
        [
            "status",
            "sequence_id",
            "bam_length",
            "fasta_length",
        ],
        mismatch_rows,
    )


def sha256_file(path):
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)

            if not chunk:
                break

            digest.update(chunk)

    return digest.hexdigest()


def write_reference_checksums(
    myo_fasta,
    myo_gtf,
    car_fasta,
    output_dir,
):
    files = [
        ("Myotis_lucifugus_genome", myo_fasta),
        ("Myotis_lucifugus_annotation", myo_gtf),
        ("Carollia_perspicillata_genome", car_fasta),
    ]

    rows = []

    for label, path in files:
        print(
            f"[Checksum] {path.name}",
            file=sys.stderr,
        )

        rows.append(
            [
                label,
                path.name,
                sha256_file(path),
                path.stat().st_size,
                str(path),
            ]
        )

    write_tsv(
        output_dir / "reference_files_sha256.tsv",
        [
            "resource",
            "filename",
            "sha256",
            "compressed_size_bytes",
            "local_source_path",
        ],
        rows,
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--myo-fa", required=True, type=Path)
    parser.add_argument("--myo-gtf", required=True, type=Path)
    parser.add_argument("--car-fa", required=True, type=Path)

    args = parser.parse_args()

    output_dir = args.repo / "provenance" / "mapping_qc"
    reference_output_dir = (
        args.repo / "provenance" / "reference_databases"
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    reference_output_dir.mkdir(parents=True, exist_ok=True)

    mapping_rows = extract_mapping_metrics(
        args.run,
        output_dir,
    )

    assigned_lookup, _ = parse_featurecounts_summary(
        args.run,
        mapping_rows,
        output_dir,
    )

    parse_count_matrix(
        args.run,
        args.myo_gtf,
        mapping_rows,
        assigned_lookup,
        output_dir,
    )

    clean_mismatch_table(
        args.run,
        output_dir,
    )

    validate_reference_dictionary(
        args.run,
        args.myo_fa,
        reference_output_dir,
    )

    write_reference_checksums(
        args.myo_fa,
        args.myo_gtf,
        args.car_fa,
        reference_output_dir,
    )

    print()
    print("[OK] Mapping-QC and reference files created.")
    print(f"[OK] Mapping QC: {output_dir}")
    print(f"[OK] References: {reference_output_dir}")


if __name__ == "__main__":
    main()
