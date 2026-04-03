#!/bin/bash
#SBATCH --job-name=igv_verify
#SBATCH --output=igv_verify_%j.log
#SBATCH --error=igv_verify_%j.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=medium

set -euo pipefail

# Activate conda env
source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
conda activate 16S_fix

# File paths
polished_assembly='/private/groups/russelllab/jodie/Jacobs_et_al_2026_de_novo_wRi_merrill_23_assembly/v3_assembly_fix_16S/manual_insertion/20260402_wRi_M23_manual_insertion_pilon.fasta'
short_r1="/private/groups/russelllab/jodie/Jacobs_et_al_2026_de_novo_wRi_merrill_23_assembly/data/short_reads/wRi_R1.fastq.gz"
short_r2="/private/groups/russelllab/jodie/Jacobs_et_al_2026_de_novo_wRi_merrill_23_assembly/data/short_reads/wRi_R2.fastq.gz"
long_reads="/private/groups/russelllab/jodie/Jacobs_et_al_2026_de_novo_wRi_merrill_23_assembly/data/basecalled/merge/20251218_Merrill23_wRi.plus_bootcamp.wRiM23.fastq.gz"

# Coordinates
insertion_start=651515
insertion_end=$((651515 + 144486))
assembly_16S_start=$((651515 + 765419 - 728866))
assembly_16S_end=$((651515 + 766923 - 728866))
padding=10000

# Output dir
outdir="/private/groups/russelllab/jodie/Jacobs_et_al_2026_de_novo_wRi_merrill_23_assembly/v3_assembly_fix_16S/manual_insertion"
cd $outdir

# Index assembly
echo "[1/6] Indexing assembly..."
samtools faidx $polished_assembly
bwa index $polished_assembly
minimap2 -d ${polished_assembly%.fasta}.mmi $polished_assembly

# Align short reads (paired)
echo "[2/6] Aligning short reads..."
bwa mem -t 16 $polished_assembly $short_r1 $short_r2 | \
    samtools sort -@ 16 -o short_reads_vs_polished.bam
samtools index short_reads_vs_polished.bam
samtools flagstat short_reads_vs_polished.bam

# Align long reads
echo "[3/6] Aligning long reads..."
minimap2 -ax map-ont \
    --secondary=no \
    -t 16 \
    $polished_assembly \
    $long_reads | \
    samtools sort -@ 16 -o long_reads_vs_polished.bam
samtools index long_reads_vs_polished.bam
samtools flagstat long_reads_vs_polished.bam

# Coverage depth files for plotting
echo "[4/6] Computing coverage depth..."
samtools depth -a short_reads_vs_polished.bam > short_read_depth.txt
samtools depth -a long_reads_vs_polished.bam  > long_read_depth.txt

# Get the actual contig name
contig_name=$(grep ">" $polished_assembly | sed 's/>//' | awk '{print $1}')
echo "Contig name: $contig_name"

# Rebuild BED with correct contig name
assembly_length=$(samtools view -H long_reads_vs_polished.bam | grep "^@SQ" | grep -oP "LN:\K[0-9]+")

cat << EOF > regions_of_interest.bed
${contig_name}	1	${assembly_length}	whole_assembly
${contig_name}	$((insertion_start - padding))	$((insertion_start + padding))	left_junction
${contig_name}	$((insertion_end - padding))	$((insertion_end + padding))	right_junction
${contig_name}	$((assembly_16S_start - padding))	$((assembly_16S_end + padding))	16S_rRNA
${contig_name}	$((insertion_start - padding))	$((insertion_end + padding))	full_insertion
EOF

cat regions_of_interest.bed

# Reindex polished assembly
samtools faidx $polished_assembly

create_report \
    regions_of_interest.bed \
    --fasta $polished_assembly \
    --tracks long_reads_vs_polished.bam short_reads_vs_polished.bam \
    --output igv_report.html

# Generate IGV report
echo "[6/6] Generating IGV report..."
create_report \
    regions_of_interest.bed \
    --fasta $polished_assembly \
    --tracks long_reads_vs_polished.bam short_reads_vs_polished.bam \
    --output igv_report.html

echo "Done! IGV report: $outdir/igv_report.html"