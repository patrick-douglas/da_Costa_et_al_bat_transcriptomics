# Comparative Host Transcriptomics in *Carollia perspicillata*

This repository contains supporting data, quality-control summaries,
reference-resource documentation, and reproducibility scripts associated with
the manuscript:

> **Comparative Host Transcriptomics Reveals Distinct Metabolic and Neural
> Programs in the Liver and Telencephalon of *Carollia perspicillata***

## Study overview

This study compares host gene-expression profiles between liver and
telencephalon tissues from six adult *Carollia perspicillata* individuals.
Each animal contributed paired liver and telencephalon samples, resulting in
12 poly(A)-enriched, single-end RNA-seq libraries.

The primary analysis focuses on comparative host transcriptomics. Viral
classification was performed as a secondary exploratory analysis of the
non-host read fraction. Viral-associated assignments are interpreted
conservatively and do not establish active infection, viral replication,
tissue colonization, or a causal relationship with host gene expression.

## Experimental design

- Six adult *Carollia perspicillata* individuals
- Six liver libraries
- Six telencephalon libraries
- Paired tissue design
- Twelve single-end RNA-seq libraries
- Poly(A)-enriched RNA libraries
- Ion GeneStudio S5 sequencing platform
- Differential-expression model accounting for animal identity and tissue

The paired differential-expression design was:

```text
~ animal + tissue
```

## Host-reference strategy

Host-read classification and gene-level alignment were performed as separate
analytical steps.

Quality-filtered reads were partitioned using
[Kraken2](https://github.com/DerrickWood/kraken2), as implemented in the
[Meta-Transcriptome Detector](https://github.com/FEI38750/MTD), against a
custom database combining:

- *Carollia perspicillata*: NCBI assembly GCA_056371365.1, mCarPer1.2
- *Myotis lucifugus*: Myoluc2.0

This combined Kraken2 database was used exclusively to classify reads as
host-associated or non-host.

Host-associated reads were subsequently aligned with
[Magic-BLAST](https://ncbi.github.io/magicblast/) against the annotated
*Myotis lucifugus* Myoluc2.0 genome from Ensembl release 115.

The exact reference files used for this stage were:

- `Myotis_lucifugus.Myoluc2.0.dna.toplevel.fa.gz`
- `Myotis_lucifugus.Myoluc2.0.115.gtf.gz`

Alignment files were converted to BAM format and coordinate-sorted with
[SAMtools](https://www.htslib.org/). Gene-level counts were generated with
[featureCounts](https://subread.sourceforge.net/) using the corresponding
Ensembl release 115 GTF annotation.

The use of *M. lucifugus* for gene-level quantification may preferentially
retain more conserved transcripts and reduce recovery of divergent or
incompletely represented *C. perspicillata* loci. Results are therefore
interpreted as comparative expression patterns recovered under this
cross-species reference framework rather than as a complete reconstruction
of the *C. perspicillata* transcriptome.

## Host mapping and annotation summary

Across the 12 libraries:

| Metric | Result |
|---|---:|
| Raw reads before fastp | 202,114,453 |
| Reads retained after fastp | 183,168,903 |
| Reads classified as host-associated | 160,211,513 |
| Pooled host-classification rate | 87.47% |
| Primary Magic-BLAST alignments | 59,419,590 |
| Primary alignment rate among host-associated reads | 37.09% |
| Primary alignment rate among quality-filtered reads | 32.44% |
| Alignments assigned to genes by featureCounts | 33,021,679 |
| Pooled featureCounts assignment rate among primary alignments | 55.57% |
| Genes represented in the Ensembl release 115 GTF | 25,849 |
| Genes represented in the featureCounts matrix | 25,849 |
| Genes detected in at least one library | 17,291 |
| Detected fraction of annotated genes | 66.89% |
| Pooled estimated base-substitution mismatch rate | 2.21% |

Across individual libraries:

| Metric | Median | Range |
|---|---:|---:|
| Host-classification rate | 84.87% | 69.93–95.43% |
| Primary alignment rate among host-associated reads | 26.86% | 16.11–64.92% |
| featureCounts assignment rate among primary alignments | 67.89% | 43.28–77.25% |
| Estimated base-substitution mismatch rate | 3.54% | 1.42–4.15% |

All primary alignments contained an `NM` tag. The estimated substitution
mismatch rate was calculated from the alignment `NM` value after subtracting
insertion and deletion lengths encoded in the CIGAR string. This metric refers
only to primary Magic-BLAST alignments against *Myotis lucifugus* Myoluc2.0.

Complete per-library values, denominators, validation results, and generation
scripts are available under
[`provenance/mapping_qc/`](provenance/mapping_qc/).

## Reference validation

The sequence dictionary extracted from a representative BAM file was compared
with the Myoluc2.0 FASTA reference used in the analysis.

The validation found:

- 11,654 reference sequences in the BAM header
- 11,654 reference sequences in the FASTA
- No sequence identifiers found exclusively in either resource
- No sequence-length mismatches
- Exact BAM-versus-FASTA dictionary correspondence

The complete validation table is available at:

[`provenance/reference_databases/reference_dictionary_validation.tsv`](provenance/reference_databases/reference_dictionary_validation.tsv)

## Repository contents

```text
.
├── FileS1.xls
├── README.md
└── provenance
    ├── README.md
    ├── mapping_qc
    │   ├── README.md
    │   ├── host_mapping_efficiency.tsv
    │   ├── host_mapping_summary.tsv
    │   ├── host_mapping_summary_by_tissue.tsv
    │   ├── featurecounts_assignment_summary.tsv
    │   ├── featurecounts_assignment_by_library.tsv
    │   ├── featurecounts_assignment_by_tissue.tsv
    │   ├── featurecounts_alignment_record_diagnostics.tsv
    │   ├── annotation_coverage_summary.tsv
    │   ├── annotation_coverage_by_library.tsv
    │   ├── annotation_coverage_descriptive_summary.tsv
    │   ├── host_substitution_mismatch_summary.tsv
    │   ├── host_substitution_mismatch_by_library.tsv
    │   └── scripts
    │       ├── build_mapping_qc_files.py
    │       └── finalize_mapping_qc_outputs.py
    └── reference_databases
        ├── reference_dictionary_validation.tsv
        ├── reference_dictionary_mismatches.tsv
        └── reference_files_sha256.tsv
```

## Data availability

Raw sequencing data are available through the NCBI Sequence Read Archive
under BioProject
[PRJNA1233355](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1233355).

The negative process-control sequencing library is available under SRA
accession
[SRR39883723](https://www.ncbi.nlm.nih.gov/sra/SRR39883723).

Large reference genomes and raw sequencing files are not duplicated in this
repository. Reference identities, filenames, validation results, and SHA-256
checksums are provided to support reconstruction and file-integrity
verification.

## Reproducibility

The scripts under
[`provenance/mapping_qc/scripts/`](provenance/mapping_qc/scripts/)
regenerate the public mapping, annotation-coverage, mismatch, and
reference-validation tables from the original analysis outputs.

The scripts require local access to:

- the original analysis output directory;
- the *M. lucifugus* Myoluc2.0 Ensembl release 115 genome FASTA;
- the corresponding Ensembl release 115 GTF;
- the *C. perspicillata* GCA_056371365.1 genome FASTA;
- Python 3;
- SAMtools.

Local absolute paths are supplied as command-line arguments and are not stored
in the public result tables.

## Supplementary file

`FileS1.xls` contains the supplementary table associated with the study.

Its description and correspondence with the final manuscript will be updated
before publication to ensure consistency with the revised figure, table, and
supplementary-file numbering.

## Citation

The manuscript is currently under revision. The complete citation and DOI
will be added after publication.

## Contact

Questions about the repository can be submitted through the GitHub
[Issues](../../issues) page.
