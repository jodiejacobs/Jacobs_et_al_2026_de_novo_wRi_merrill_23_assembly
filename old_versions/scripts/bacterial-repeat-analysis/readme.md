# Bacterial Repeat Analysis Pipeline

This Snakemake pipeline analyzes bacterial genomes for repetitive elements, mobile elements, insertion sequences, and other repeat structures. It's designed to work with bacterial genomes like the Wolbachia or other bacterial symbionts.

## Features

- **RepeatMasker analysis** to identify known repetitive elements
- **ISEScan** for identification of bacterial insertion sequences (IS)
- **RED (Repeat Explorer)** for de novo repeat identification
- **Custom interspersed repeat detection** using self-BLAST
- **Visualization tools** for repeat distribution in bacterial genomes
- **Comprehensive HTML report** with summary statistics and figures

## Requirements

- Linux environment with Conda/Mamba
- Snakemake
- Input bacterial genome assembly in FASTA format

## Directory Structure

```
├── config/
│   └── config.yaml           # Pipeline configuration file
├── envs/                     # Conda environment files
│   ├── repeatmasker.yaml
│   ├── isescan.yaml
│   ├── red.yaml
│   ├── blast.yaml
│   ├── report.yaml
│   └── plotting.yaml
├── scripts/                  # Analysis scripts
│   ├── filter_self_blast.py
│   ├── generate_report.py
│   └── plot_genome.py
├── Snakefile                 # Main Snakemake workflow
├── README.md                 # This file
└── data/                     # Input data (not included in repo)
    └── genomes/
        └── your-genome.fasta
```

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/bacterial-repeat-analysis.git
cd bacterial-repeat-analysis
```

2. Create a base environment with Snakemake:
```bash
mamba create -n snakemake -c conda-forge -c bioconda snakemake-minimal
conda activate snakemake
```

3. The pipeline will automatically create tool-specific environments as needed.

## Usage

1. Place your bacterial genome assembly in the `data/genomes/` directory.

2. Edit the configuration file `config/config.yaml` to set your input genome path, output directory, and analysis parameters.

3. Run the pipeline:
```bash
snakemake --use-conda -j <number_of_cores>
```

For a dry run to check the workflow:
```bash
snakemake -n --use-conda
```

### SLURM Cluster Execution

To run on a SLURM cluster:

```bash
snakemake --use-conda --profile slurm --cluster-config cluster.yaml
```

A sample `cluster.yaml` is included in the repository.

## Output Files

The pipeline generates the following outputs in the specified output directory:

### RepeatMasker Output
- `repeatmasker/your-genome.masked`: Soft-masked genome sequence
- `repeatmasker/your-genome.out`: Detailed list of identified repeats
- `repeatmasker/your-genome.tbl`: Summary table of repeat categories

### ISEScan Output
- `isescan/output/your-genome.features`: Detailed list of insertion sequences
- Other ISEScan annotation files

### RED Output
- `red/red_repeats.txt`: De novo identified repeats

### Interspersed Repeat Analysis
- `interspersed_repeats/repeats.gff`: GFF file with interspersed repeats

### Visualization and Reports
- `plots/genome_with_repeats.pdf`: Circular genome plot with repeat elements
- `plots/genome_with_repeats_linear.pdf`: Linear genome plot
- `plots/genome_with_repeats_density.pdf`: Repeat density plot
- `repeat_summary_report.html`: Comprehensive HTML report

## Customization

You can customize the analysis by modifying the `config.yaml` file. The main parameters include:

- `assembly`: Path to your genome assembly
- `output_dir`: Where results will be stored
- `threads`: Number of CPU threads to use
- `min_repeat_length`: Minimum length for interspersed repeats
- `min_repeat_identity`: Minimum sequence identity for repeats

## Troubleshooting

If you encounter any issues:

1. Check the log files in the `logs/` directory
2. Make sure all dependencies are installed
3. For RepeatMasker, ensure the bacterial RepBase is correctly configured

## Citation

If you use this pipeline in your research, please cite the following tools:

- **RepeatMasker**: Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-4.0. 2013-2015
- **ISEScan**: Xie Z, Tang H. ISEScan: automated identification of insertion sequence elements in prokaryotic genomes. Bioinformatics. 2017
- **RED**: Girgis HZ. Red: an intelligent, rapid, accurate tool for detecting repeats de-novo on the genomic scale. BMC Bioinformatics. 2015

## License

This pipeline is available under the MIT License. Feel free to modify and distribute it.
