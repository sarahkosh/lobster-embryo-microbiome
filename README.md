# American Lobster Embryo and Larval Microbiome Analysis

R scripts for processing and analyzing 16S rRNA amplicon sequencing data from American lobster (*Homarus americanus*) embryos, recently hatched larvae, and tank water.

## Associated Publication

Koshak, J. S., Song, B., Jellison, B., et al. (2026). The bacterial community composition of American lobster (*Homarus americanus*) embryos and recently hatched larvae held under experimental laboratory conditions. *Scientific Reports*, 16, 9045.

https://doi.org/10.1038/s41598-026-35387-8

## Data Availability

Raw sequencing data are publicly available through the NCBI Sequence Read Archive (SRA) under BioProject accession **PRJNA1227997**.

[View the BioProject on NCBI](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1227997)

## Repository Contents

| File | Description |
|---|---|
| `AdditionalFile1.R` | DADA2 workflow for processing 16S rRNA amplicon sequencing data from lobster embryos. |
| `AdditionalFile2.R` | DADA2 workflow for processing 16S rRNA amplicon sequencing data from lobster larvae and tank water samples. |
| `AdditionalFile3.R` | Taxonomic analysis and visualization using the `phyloseq` R package. |
| `AdditionalFile4.R` | Core microbiome analysis. |

## Taxonomic Classification

Taxonomic assignments were performed using the SILVA rRNA database, version 132.

## Reproducibility

These scripts document the bioinformatics and downstream analytical workflows used in the associated publication. Refer to the publication for details on the experimental design, sequencing methods, and statistical analyses.
