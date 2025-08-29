#!/bin/bash

set -e

if [ "$#" -ne 6 ]; then
    echo "Usage: bash circularize_assembly.sh assembly.fasta bridge_read.fasta contig_8 contig_18 output.fasta"
    exit 1
fi

assembly_fasta=$1
bridge_fasta=$2
contig_8=$3
contig_18=$4
# contig_36=$5
output_fasta=$6

read_id=$(grep ">" "$bridge_fasta" | head -n 1 | sed 's/>//')

echo "🔄 Aligning bridge read to assembly..."
minimap2 -ax map-ont "$assembly_fasta" "$bridge_fasta" > bridge_read.sam
samtools view -bS bridge_read.sam | samtools sort -o bridge_read.sorted.bam
samtools index bridge_read.sorted.bam

echo "📌 Identifying read positions..."
positions=$(samtools view bridge_read.sorted.bam | grep "$read_id" | awk '{print $3, $4, length($10)}')

declare -A contig_positions
while read -r contig start read_length; do
    contig_positions["$contig"]="$start $((start + read_length))"
done <<< "$positions"

for contig in "$contig_8" "$contig_18" "$contig_36"; do
    if [[ -z "${contig_positions[$contig]}" ]]; then
        echo "❌ Error: The read does not map to $contig."
        exit 1
    fi
done

contig_8_start=$(echo "${contig_positions[$contig_8]}" | cut -d' ' -f1)
contig_18_end=$(echo "${contig_positions[$contig_18]}" | cut -d' ' -f2)
# contig_36_end=$(echo "${contig_positions[$contig_36]}" | cut -d' ' -f2)

contig_8_len=$(samtools faidx "$assembly_fasta" | grep -P "^$contig_8\s" | cut -f2)
contig_18_len=$(samtools faidx "$assembly_fasta" | grep -P "^$contig_18\s" | cut -f2)
contig_36_len=$(samtools faidx "$assembly_fasta" | grep -P "^$contig_36\s" | cut -f2)

echo "✂️ Extracting relevant contig regions..."
samtools faidx "$assembly_fasta" "$contig_8":1-"$contig_8_start" > trimmed_contig_8.fasta
samtools faidx "$assembly_fasta" "$contig_18":"$contig_18_end"-"$contig_18_len" > trimmed_contig_18.fasta
# samtools faidx "$assembly_fasta" "$contig_36":"$contig_36_end"-"$contig_36_len" > trimmed_contig_36.fasta

echo "✂️ Extracting non-overlapping bridge read region..."
read_len=$(samtools faidx "$bridge_fasta" "$read_id" | awk 'NR==2 {print length($0)}')

if [[ $contig_8_end -lt $read_len && $contig_18_start -gt 1 ]]; then
    samtools faidx "$bridge_fasta" "$read_id":"$contig_8_end"-"$contig_18_start" > trimmed_bridge_read.fasta
else
    echo "❌ Error: Invalid bridge read extraction coordinates"
    exit 1
fi

echo "🔗 Merging sequences to form circularized genome..."
cat trimmed_contig_8.fasta trimmed_bridge_read.fasta trimmed_contig_36.fasta trimmed_contig_18.fasta > circularized.fasta

echo "✨ Polishing the circularized genome..."
medaka_consensus -i "$bridge_fasta" -d circularized.fasta -o polished_output -m r941_min_hac_g507

cp polished_output/consensus.fasta "$output_fasta"
echo "🎉 Circularized assembly saved as $output_fasta"