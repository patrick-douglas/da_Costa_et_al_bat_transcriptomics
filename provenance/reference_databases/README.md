# Host reference resources

This directory documents the host genomic resources used during read
classification, alignment, and gene-level quantification.

## *Carollia perspicillata*

- Assembly accession: GCA_056371365.1
- Assembly name: mCarPer1.2
- Genome file:
  `GCA_056371365.1_mCarPer1.2_genomic.fna.gz`
- Role in the analysis: included in the combined Kraken2 database used for
  host-read classification.

The *C. perspicillata* genome was not used as the final gene-level alignment
and quantification reference because the analysis required a matched
genome-annotation resource compatible with the downstream workflow.

## *Myotis lucifugus*

- Assembly: Myoluc2.0
- Source: Ensembl release 115
- Genome file:
  `Myotis_lucifugus.Myoluc2.0.dna.toplevel.fa.gz`
- Annotation file:
  `Myotis_lucifugus.Myoluc2.0.115.gtf.gz`

The Myoluc2.0 genome was used in two contexts:

1. as one component of the combined Kraken2 host-classification database;
2. as the annotated reference used for Magic-BLAST alignment and
   featureCounts gene-level quantification.

## Validation files

### `reference_dictionary_validation.tsv`

Compares sequence identifiers and sequence lengths in a representative BAM
header against the Myoluc2.0 FASTA sequence dictionary.

The validation found:

- 11,654 reference sequences in the BAM;
- 11,654 reference sequences in the FASTA;
- no sequence identifiers exclusive to either resource;
- no sequence-length mismatches;
- an exact dictionary match.

### `reference_dictionary_mismatches.tsv`

Lists sequence-identifier or sequence-length mismatches, when present. The
file contains only its header when no mismatch is detected.

### `reference_files_sha256.tsv`

Reports SHA-256 checksums and compressed file sizes for the principal
*Carollia perspicillata* and *Myotis lucifugus* reference resources.

The genome and annotation files themselves are not redistributed here because
of their size. The checksums permit verification of locally obtained copies.
