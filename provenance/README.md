# Analysis provenance

This directory contains quality-control results, reference-resource
documentation, checksums, and reproducibility scripts supporting the
comparative host-transcriptome analysis.

## Contents

### [`mapping_qc/`](mapping_qc/)

Contains per-library, tissue-level, and pooled summaries for:

- read retention after fastp processing;
- host-read classification using the combined *Carollia perspicillata* and
  *Myotis lucifugus* Kraken2 database;
- primary Magic-BLAST alignment against *Myotis lucifugus* Myoluc2.0;
- featureCounts assignment to annotated genes;
- annotation coverage;
- estimated base-substitution mismatch rates.

The directory also contains the scripts used to generate and finalize the
public quality-control tables.

### [`reference_databases/`](reference_databases/)

Contains:

- reference filenames and SHA-256 checksums;
- validation of the BAM sequence dictionary against the
  *Myotis lucifugus* Myoluc2.0 FASTA;
- documentation of the role of each host reference resource.

Large genome, annotation, alignment, and sequencing files are not
redistributed in this repository. Their identities, checksums, and validation
results are provided to support reconstruction and file-integrity
verification.
