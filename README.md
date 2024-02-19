# GenoDiplo: Advanced Genomics and Comparative Genomics Analysis Pipelines for Diplomonads
![GenoDiplo Logo](logo.png "GenoDiplo Logo")

## Overview
GenoDiplo offers an integrated workflow for the assembly and annotation of diplomonad genomes, leveraging hybrid sequencing technologies. This repository encompasses tools and protocols for preprocessing, assembly, structural and functional annotation, and comparative genomics analysis.

### Hybrid Genome Assembly Pipeline
The Hybrid Genome and Transcriptome Assembly workflow integrates pre-assembly and post-assembly stages, utilizing a hybrid approach that combines long and short DNA reads. The process includes:

- **Preprocessing:** Quality checks, read statistics, and contaminant removal.
- **Assembly:** Utilizes Masurca and FLYE for constructing the genome, followed by Pilon for assembly polishing.
- **Evaluation:** Validation of the transcriptome, quality assessments, and coverage calculations.

### Structural and Functional Annotation Pipeline
This pipeline addresses the unique challenges of diplomonad genomes, such as the scarcity of introns, by combining Prodigal and GlimmerHMM for annotation:

- **Gene Prediction:** Filtering overlapping gene predictions between Prodigal and GlimmerHMM.
- **Functional Annotation:** Utilizes BLASTp, eggNOG-mapper, and InterProScan for comprehensive annotation.

### Comparative Genomics Analysis

#### Genome Structure Level
- **Repetitive Region Identification:** Employing RepeatMasker with a custom repeat library created by RepeatModeler.
- **tRNA Annotation:** Using tRNAscan-SE and Stainglass for tRNA gene and tandem repeat visualization.

#### Protein Domain Level
- **Evolutionary Analysis:** Integrating IPR superfamily entries with OrthoFinder to explore evolutionary relationships and functional characteristics.
- **Superfamily Heatmaps:** Visualizing expanded superfamilies to highlight adaptive and evolutionary pressures.

#### Sequence Similarity Level
- **Orthologous Groups Analysis:** Using UpSet plots to illustrate overlaps among species, emphasizing shared genetic content and horizontal gene transfer.

## Getting Started
Please refer to the [Installation](#installation) and [Usage](#usage) sections for details on how to deploy and use GenoDiplo in your projects.

