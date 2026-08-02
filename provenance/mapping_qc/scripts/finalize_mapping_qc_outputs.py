#!/usr/bin/env python3

import csv
import statistics
from pathlib import Path


REPO = Path(__file__).resolve().parents[3]
MAPPING_DIR = REPO / "provenance" / "mapping_qc"
REFERENCE_DIR = REPO / "provenance" / "reference_databases"


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writerow(header)
        writer.writerows(rows)


def percent(numerator, denominator):
    if denominator == 0:
        return "NA"

    return f"{100.0 * numerator / denominator:.4f}"


def finalize_featurecounts():
    original = MAPPING_DIR / "featurecounts_assignment_summary.tsv"
    diagnostics = (
        MAPPING_DIR / "featurecounts_alignment_record_diagnostics.tsv"
    )

    if diagnostics.exists():
        rows = read_tsv(diagnostics)
    else:
        rows = read_tsv(original)

        required = {
            "sample",
            "featurecounts_input_records",
            "featurecounts_assigned_reads",
            "primary_mapped_reads",
        }

        if not rows or not required.issubset(rows[0]):
            raise RuntimeError(
                "The original featureCounts diagnostic table "
                "could not be identified."
            )

        write_tsv(
            diagnostics,
            list(rows[0].keys()),
            [
                [row[column] for column in rows[0].keys()]
                for row in rows
            ],
        )

    clean_rows = []
    percentages = []

    for row in rows:
        sample = row["sample"]
        primary = int(row["primary_mapped_reads"])
        assigned = int(row["featurecounts_assigned_reads"])
        assignment_percent = 100.0 * assigned / primary

        tissue = (
            "Liver"
            if "_LIVER_" in sample
            else "Telencephalon"
            if "_TEL_" in sample
            else "Unknown"
        )

        percentages.append(assignment_percent)

        clean_rows.append(
            {
                "sample": sample,
                "tissue": tissue,
                "primary_mapped_alignments": primary,
                "featurecounts_assigned_alignments": assigned,
                "assigned_percent_of_primary_mapped_alignments": (
                    assignment_percent
                ),
            }
        )

    write_tsv(
        MAPPING_DIR / "featurecounts_assignment_by_library.tsv",
        [
            "sample",
            "tissue",
            "primary_mapped_alignments",
            "featurecounts_assigned_alignments",
            "assigned_percent_of_primary_mapped_alignments",
        ],
        [
            [
                row["sample"],
                row["tissue"],
                row["primary_mapped_alignments"],
                row["featurecounts_assigned_alignments"],
                (
                    f'{row["assigned_percent_of_primary_mapped_alignments"]:.4f}'
                ),
            ]
            for row in clean_rows
        ],
    )

    total_primary = sum(
        row["primary_mapped_alignments"]
        for row in clean_rows
    )

    total_assigned = sum(
        row["featurecounts_assigned_alignments"]
        for row in clean_rows
    )

    summary_rows = [
        ["libraries", len(clean_rows)],
        ["primary_mapped_alignments", total_primary],
        ["featurecounts_assigned_alignments", total_assigned],
        [
            "pooled_assignment_percent_of_primary_mapped_alignments",
            percent(total_assigned, total_primary),
        ],
        [
            "assignment_percent_mean",
            f"{statistics.mean(percentages):.4f}",
        ],
        [
            "assignment_percent_median",
            f"{statistics.median(percentages):.4f}",
        ],
        [
            "assignment_percent_minimum",
            f"{min(percentages):.4f}",
        ],
        [
            "assignment_percent_maximum",
            f"{max(percentages):.4f}",
        ],
    ]

    write_tsv(
        original,
        ["metric", "value"],
        summary_rows,
    )

    tissue_rows = []

    for tissue in ("Liver", "Telencephalon"):
        subset = [
            row
            for row in clean_rows
            if row["tissue"] == tissue
        ]

        values = [
            row["assigned_percent_of_primary_mapped_alignments"]
            for row in subset
        ]

        primary = sum(
            row["primary_mapped_alignments"]
            for row in subset
        )

        assigned = sum(
            row["featurecounts_assigned_alignments"]
            for row in subset
        )

        tissue_rows.append(
            [
                tissue,
                len(subset),
                primary,
                assigned,
                percent(assigned, primary),
                f"{statistics.mean(values):.4f}",
                f"{statistics.median(values):.4f}",
                f"{min(values):.4f}",
                f"{max(values):.4f}",
            ]
        )

    write_tsv(
        MAPPING_DIR / "featurecounts_assignment_by_tissue.tsv",
        [
            "tissue",
            "libraries",
            "primary_mapped_alignments",
            "featurecounts_assigned_alignments",
            "pooled_assignment_percent",
            "mean_assignment_percent",
            "median_assignment_percent",
            "minimum_assignment_percent",
            "maximum_assignment_percent",
        ],
        tissue_rows,
    )


def summarize_annotation_coverage():
    path = MAPPING_DIR / "annotation_coverage_by_library.tsv"
    rows = read_tsv(path)

    values = [
        float(row["detected_gene_percent_of_annotated_genes"])
        for row in rows
    ]

    write_tsv(
        MAPPING_DIR / "annotation_coverage_descriptive_summary.tsv",
        ["metric", "value"],
        [
            ["libraries", len(rows)],
            [
                "detected_gene_percent_mean",
                f"{statistics.mean(values):.4f}",
            ],
            [
                "detected_gene_percent_median",
                f"{statistics.median(values):.4f}",
            ],
            [
                "detected_gene_percent_minimum",
                f"{min(values):.4f}",
            ],
            [
                "detected_gene_percent_maximum",
                f"{max(values):.4f}",
            ],
        ],
    )


def sanitize_metric_file(path, metrics):
    rows = read_tsv(path)

    for row in rows:
        if row.get("metric") in metrics:
            row["value"] = Path(row["value"]).name

    write_tsv(
        path,
        ["metric", "value"],
        [
            [row["metric"], row["value"]]
            for row in rows
        ],
    )


def sanitize_reference_checksums():
    path = REFERENCE_DIR / "reference_files_sha256.tsv"
    rows = read_tsv(path)

    write_tsv(
        path,
        [
            "resource",
            "filename",
            "sha256",
            "compressed_size_bytes",
        ],
        [
            [
                row["resource"],
                row["filename"],
                row["sha256"],
                row["compressed_size_bytes"],
            ]
            for row in rows
        ],
    )


def main():
    finalize_featurecounts()
    summarize_annotation_coverage()

    sanitize_metric_file(
        MAPPING_DIR / "annotation_coverage_summary.tsv",
        {
            "gtf_file",
            "featurecounts_matrix",
        },
    )

    sanitize_metric_file(
        REFERENCE_DIR / "reference_dictionary_validation.tsv",
        {
            "bam_used",
            "fasta_used",
        },
    )

    sanitize_reference_checksums()

    log_file = MAPPING_DIR / "build_mapping_qc_files.log"

    if log_file.exists():
        log_file.unlink()

    print("[OK] Final public mapping-QC files prepared.")
    print(f"[OK] Repository: {REPO}")


if __name__ == "__main__":
    main()
