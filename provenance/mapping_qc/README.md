# Host mapping and annotation quality control

This directory contains per-library and pooled quality-control statistics for
the comparative host-transcriptome analysis.

## Reference strategy

Host-read classification and gene-level alignment were performed as separate
analytical steps.

Quality-filtered reads were first classified using a custom Kraken2 database
combining genomic sequences from:

- *Carollia perspicillata*: NCBI assembly GCA_056371365.1, mCarPer1.2
- *Myotis lucifugus*: Myoluc2.0

The combined Kraken2 database was used exclusively to partition reads into
host-associated and non-host fractions.

Host-associated reads were subsequently aligned with Magic-BLAST against the
annotated *Myotis lucifugus* Myoluc2.0 genome. Gene-level quantification was
performed with featureCounts using the corresponding Ensembl release 115 GTF
annotation.

## Files

### `host_mapping_efficiency.tsv`

Per-library read counts and percentages for:

- reads before fastp processing;
- reads retained after fastp;
- reads classified as host-associated by the combined Kraken2 database;
- primary Magic-BLAST alignments against *M. lucifugus*;
- host-classification and alignment percentages.

### `host_mapping_summary.tsv`

Pooled and library-level descriptive statistics.

### `host_mapping_summary_by_tissue.tsv`

Descriptive statistics calculated separately for liver and telencephalon
libraries.

### `featurecounts_assignment_summary.tsv`

Per-library featureCounts assignment statistics.

The field `assigned_percent_of_primary_mapped_reads` represents the fraction
of primary Magic-BLAST alignments assigned to annotated gene features.

### `annotation_coverage_by_library.tsv`

Per-library numbers and percentages of annotated genes with nonzero counts,
together with featureCounts assignment rates.

### `annotation_coverage_summary.tsv`

Overall number of genes represented in the Ensembl release 115 GTF,
genes present in the featureCounts matrix, and genes detected in at least one
library.

### `host_substitution_mismatch_by_library.tsv`

Estimated base-substitution mismatch rates among primary Magic-BLAST
alignments.

Substitution mismatches were estimated from the NM tag after subtracting
inserted and deleted bases encoded in each CIGAR string. Long-gap edit-distance
percentages were not included because they can be dominated by deletion or
gap operations and should not be interpreted as simple nucleotide mismatch
rates.

### `host_substitution_mismatch_summary.tsv`

Library-level descriptive statistics and pooled base-weighted substitution
mismatch rate.

### `scripts/build_mapping_qc_files.py`

Reproducible script used to generate the files in this directory.

## Denominators

The host-classification percentage uses quality-filtered reads as the
denominator.

The primary-alignment percentage is reported using both:

1. host-associated reads as the denominator; and
2. all quality-filtered reads as the denominator.

Mismatch statistics refer only to primary Magic-BLAST alignments against
*Myotis lucifugus* Myoluc2.0.
