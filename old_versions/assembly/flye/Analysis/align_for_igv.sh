# Map all reads to genome assembly for visualization in IGV 
minimap2 -ax map-ont --secondary=no --cs --eqx -t 16 \
/private/groups/russelllab/jodie/merrill_23_wRi_genome/flye/merge_bootcamp_data/assembly.fasta \
/private/groups/russelllab/jodie/merrill_23_wRi_genome/flye/merge_bootcamp_data/merrill_23_wRi_filtered.fastq.gz | \
samtools view -b -q 20 | \
samtools sort -o /private/groups/russelllab/jodie/merrill_23_wRi_genome/flye/merge_bootcamp_data/Analysis/aligned_reads.bam
samtools index /private/groups/russelllab/jodie/merrill_23_wRi_genome/flye/merge_bootcamp_data/Analysis/aligned_reads.bam

## Get spanning reads
samtools view -H /private/groups/russelllab/jodie/merrill_23_wRi_genome/flye/merge_bootcamp_data/Analysis/aligned_reads.bam > header.sam
samtools view /private/groups/russelllab/jodie/merrill_23_wRi_genome/flye/merge_bootcamp_data/Analysis/aligned_reads.bam | grep "SA:" > body.sam
cat header.sam body.sam | samtools view -Sb - > /private/groups/russelllab/jodie/merrill_23_wRi_genome/flye/merge_bootcamp_data/Analysis/spanning_reads.bam
samtools index /private/groups/russelllab/jodie/merrill_23_wRi_genome/flye/merge_bootcamp_data/Analysis/spanning_reads.bam

# Extract reads that map to both contig 8 and 18 
bedtools bamtobed -i aligned_reads.bam > aligned_reads.bed
# Extract read names
cut -f 1,4 aligned_reads.bed > read_contig_pairs.txt
#Identify reads appearing in mulitple contigs
sort -k2,2 read_contig_pairs.txt | \
uniq -c | \
awk '$2 > 1 {print $3}' > reads_in_multiple_contigs.txt


# Extract the reads in multiple contigs
grep -Ff reads_in_multiple_contigs.txt aligned_reads.bed > reads_spanning_multiple_contigs.bed

samtools view -N reads_in_multiple_contigs.txt -h aligned_reads.bam > reads_spanning_multiple_contigs.sam
samtools view -Sb reads_spanning_multiple_contigs.sam > reads_spanning_multiple_contigs.bam
samtools index reads_spanning_multiple_contigs.bam

