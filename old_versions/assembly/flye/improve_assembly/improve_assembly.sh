# Remove contigs less than 1000 bp
seqkit seq -m 1000 -g /private/groups/russelllab/jodie/merrill_23_wRi_genome/flye/merge_bootcamp_data/assembly.fasta > /private/groups/russelllab/jodie/merrill_23_wRi_genome/flye/merge_bootcamp_data/filtered_assembly.fasta
# Eliminate contigs that are from contamination 
# kraken2 --db kraken_db --threads 16 --output kraken_output.txt --report kraken_report.txt /private/groups/russelllab/jodie/merrill_23_wRi_genome/flye/merge_bootcamp_data/filtered_assembly.fasta
blastn -query filtered_assembly.fasta -db nt -out contig_species_blast.txt -outfmt "6 qseqid staxids sscinames pident length evalue bitscore" -max_target_seqs 1 -num_threads 16
awk '{print $1 "\t" $3}' contig_species_blast.txt > contig_species_mapping.txt

#seqkit grep -v -f contaminant_ids.txt assembly.fasta > cleaned_assembly.fasta

# minimap2 -ax map-ont -t 16 assembly.fasta reads.fastq | samtools sort -o aligned_reads.bam
# samtools index aligned_reads.bam

# racon -t 16 reads.fastq aligned_reads.bam assembly.fasta > polished_assembly.fasta

# medaka_consensus -i reads.fastq -d assembly.fasta -o polished_assembly.fasta -t 16

# busco -i polished_assembly.fasta -l your_lineage_odb10 -m genome -o busco_output