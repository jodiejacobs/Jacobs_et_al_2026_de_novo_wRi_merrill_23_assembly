# Complete de novo assembly of Wolbachia endosymbiont of contemporary Drosophila simulans using long-read genome sequencing

[![DOI](https://img.shields.io/badge/DOI-PRJNA1312834-blue)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1312834)

This repository contains the analysis pipeline and code for the complete de novo assembly of Wolbachia pipientis strain wRi Merrill 23, an alphaproteobacterial endosymbiont of *Drosophila simulans* collected at UC Santa Cruz in 2023.

## Background

Wolbachia pipientis infects diverse arthropods and nematodes, manipulating host phenotypes through cytoplasmic incompatibility (CI), male killing, and fertility rescue. The Riverside strain (wRi) was first identified in California *Drosophila simulans* in the 1980s and rapidly spread statewide due to exceptionally strong CI. Despite its significance in shaping *D. simulans* populations, a modern wRi genome had not been assembled, with the existing reference genome reflecting wRi from 1984.

This project provides the first contemporary complete de novo assembly of wRi from *D. simulans* collected in 2023, offering an updated genomic reference for future evolutionary and functional studies.

## Key Results

- **Complete circular assembly**: 1.26 Mb genome with 30X coverage
- **High quality**: 99.2% BUSCO completeness score
- **Modern reference**: First assembly from contemporary wRi (2023 vs 1984)
- **Comprehensive annotation**: 1,283 genes including 1,246 CDSs and 37 RNA genes

## Repository Structure

```
├── README.md
├── Snakefile                    # Main Snakemake workflow
├── basecaller-slurm.sh         # SLURM script for nanopore basecalling
├── config/
│   └── read_filtering.yaml     # Conda environment configuration
├── data/
│   ├── pod5/                   # Raw nanopore sequencing data
│   ├── aligned/                # Aligned reads
│   ├── basecalled/             # Basecalled fastq files
│   ├── flye/                   # Flye assembly outputs
│   ├── polished/               # Final polished assemblies
│   ├── short_reads/            # Illumina short reads for polishing
│   └── reference/              # D. simulans reference genome
└── busco/                      # BUSCO quality assessment results
```

## Methods Overview

### Sequencing and Data Collection
- **Sample**: wRi-infected *D. simulans* embryos from UCSC Merrill gardens (October 2023)
- **DNA extraction**: Wizard HMW DNA Extraction kit (Promega)
- **Library prep**: Native Barcoding Kit V14 (Oxford Nanopore)
- **Sequencing**: MinION Mk1B with R10 flow cell, 20 hours runtime
- **Reads generated**: 5.4M reads with adaptive host depletion

### Assembly Pipeline
1. **Basecalling**: Oxford Nanopore Dorado v0.7.3 with HAC model
2. **Read filtering**: Exclude host sequences, retain reads >3kb
3. **Assembly**: Flye assembler optimized for bacterial genomes
4. **Polishing**: POLCA with Illumina short reads
5. **Quality assessment**: BUSCO with rickettsiales_odb10 database
6. **Annotation**: Prokka v1.1.1 for gene prediction

## Requirements

### Software Dependencies
- **Nanopore tools**: Dorado v0.7.3, MinKNOW v23.07.8
- **Assembly**: Flye, POLCA, Pilon
- **Read processing**: samtools, sambamba, seqkit, BWA
- **Quality control**: BUSCO v5.7.0
- **Annotation**: Prokka v1.1.1
- **Visualization**: Proksee v1.1.2
- **Workflow management**: Snakemake

### System Requirements
- Linux-based HPC environment with SLURM scheduler
- GPU support for Dorado basecalling (4 GPUs recommended)
- Minimum 100GB RAM for assembly steps
- Python 3.8+ with mamba/conda for environment management

## Usage

### 1. Environment Setup

Create the conda environment:
```bash
mamba env create -f config/read_filtering.yaml
```

Additional environments needed:
```bash
# For assembly pipeline
mamba create -n assembly flye polca pilon bwa samtools sambamba seqkit

# For quality assessment  
mamba create -n busco busco

# For Snakemake
mamba create -n snakemake snakemake
```

### 2. Data Preparation

Place your data in the appropriate directories:
- Raw nanopore data (POD5 format) in `data/pod5/`
- Illumina short reads in `data/short_reads/`
- *D. simulans* reference genome in `data/reference/`

### 3. Running the Pipeline

Navigate to the project directory and activate the Snakemake environment:
```bash
cd /path/to/Jacobs_et_al_2026_de_novo_wRi_merrill_23_assembly/
mamba activate snakemake
```

Run the complete pipeline:
```bash
snakemake --executor slurm \
          --default-resources slurm_partition=medium runtime=720 mem_mb=1000000 \
          -j 10 \
          -s Snakefile
```

For manual basecalling (alternative approach):
```bash
sbatch basecaller-slurm.sh data/pod5/ data/aligned/ data/reference/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna.gz
```

### 4. Outputs

The final polished assembly will be available at:
```
data/polished/{sample}_wRi_M23.assembly.fasta
```

Quality assessment results will be in:
```
busco/{sample}/wRi/polished/short_summary.specific.rickettsiales_odb10.txt
```

## Key Features

- **Adaptive sequencing**: Real-time host depletion during nanopore sequencing
- **Hybrid approach**: Long-read assembly with short-read polishing
- **Quality control**: Comprehensive BUSCO assessment
- **Reproducible workflow**: Complete Snakemake pipeline with SLURM integration
- **Modular design**: Separate rules for each major processing step

## Assembly Statistics

| Metric | Value |
|--------|-------|
| Assembly length | 1,259,726 bp |
| GC content | 35.22% |
| Total genes | 1,283 |
| Protein-coding genes | 1,246 |
| RNA genes | 37 |
| rRNA genes | 3 (5S, 16S, 23S) |
| tRNA genes | 34 |
| BUSCO completeness | 99.2% |

## Data Availability

- **Raw sequencing data**: NCBI SRA under BioProject [PRJNA1312834](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1312834)
- **Assembled genome**: GenBank under BioProject PRJNA1312834
- **Analysis code**: This GitHub repository

## Citation

If you use this assembly or pipeline in your research, please cite:

```
Jacobs J, Lum A, Lee D, Gutierrez E, Dionisio J, Morey C, Mirchandani C, Sylvester L, 
Nakamoto A, Loucks H, Wanke C, Cisneros A, Shanks C, Calicchio A, Enstrom A, 
Michelon C, Okamoto F, Heath H, Malukhina K, Russell P, Nag S, Gillespie T, 
Sobolewski W, Truong Z, Russell SL. Complete de novo assembly of Wolbachia 
endosymbiont of contemporary Drosophila simulans using long-read genome sequencing. 
(preprint). 2025.
```

## Contributing

This pipeline was developed as part of ongoing research in the Russell Lab at UC Santa Cruz. For questions about the methods or to report issues, please open an issue on this repository.

## Acknowledgments

The authors acknowledge the University of California Santa Cruz Genomics Institute for providing computational resources and support for this project. The authors thank James Letcinger for the use of his computer for nanopore sequencing. Funding for this project was provided by NIH (T32 HG012344) awarded to JJ, CW, HL, AN, CS, and AS, NIH awards to SLR (R00GM135583, R35GM157189), and the NSF-GRFP awarded to AN.

## Contact

- **Corresponding Author**: Shelbi L Russell (shelbilrussell@gmail.com)
- **Lead Author**: Jodie Jacobs (jacobs.jodiem@gmail.com)
- **Lab Website**: [Russell Lab at UC Santa Cruz](https://russelllab.sites.ucsc.edu/)
