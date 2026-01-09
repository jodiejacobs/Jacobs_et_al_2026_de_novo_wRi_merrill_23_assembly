#!/bin/bash
#SBATCH --job-name=collect_assembly_results
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --dependency=singleton
#SBATCH --output=collect_results.log
#SBATCH --partition=medium

WORKDIR="multi_assembly_comparison"
cd ${WORKDIR}

MIN_SIZE=1200000
MAX_SIZE=1400000

# Initialize summary file
echo "Assembly Comparison for wRi Merrill 23" > assembly_summary.txt
echo "Date: $(date)" >> assembly_summary.txt
echo "Expected chromosome size: ${MIN_SIZE}-${MAX_SIZE} bp" >> assembly_summary.txt
echo "======================================" >> assembly_summary.txt
echo "" >> assembly_summary.txt

# Collect individual summaries
for summary in *_summary.txt; do
    if [ -f ${summary} ]; then
        cat ${summary} >> assembly_summary.txt
        echo "" >> assembly_summary.txt
    fi
done

# Create comparison table
echo "" >> assembly_summary.txt
echo "=== SUMMARY TABLE ===" >> assembly_summary.txt
echo "Assembler | N_Contigs | Total_Size | N50 | Candidate_Chromosomes" >> assembly_summary.txt
echo "----------|-----------|------------|-----|----------------------" >> assembly_summary.txt

for asm_dir in flye_default flye_conservative canu_assembly raven_assembly miniasm_assembly hifiasm_assembly; do
    if [ -d ${asm_dir} ]; then
        # Find assembly file
        if [ -f ${asm_dir}/assembly.fasta ]; then
            asm_file=${asm_dir}/assembly.fasta
        elif [ -f ${asm_dir}/wRi.contigs.fasta ]; then
            asm_file=${asm_dir}/wRi.contigs.fasta
        else
            continue
        fi
        
        # Get stats
        n_contigs=$(grep -c ">" ${asm_file} || echo "0")
        
        if [ ${n_contigs} -gt 0 ]; then
            stats=$(seqkit stats -T ${asm_file} | tail -n 1)
            total=$(echo ${stats} | awk '{print $5}')
            n50=$(echo ${stats} | awk '{print $14}')
            
            # Count candidates
            candidates=$(seqkit fx2tab -l -n ${asm_file} | awk -v min=${MIN_SIZE} -v max=${MAX_SIZE} '$2 >= min && $2 <= max' | wc -l)
            
            echo "${asm_dir} | ${n_contigs} | ${total} | ${n50} | ${candidates}" >> assembly_summary.txt
        fi
    fi
done

echo ""
echo "======================================"
echo "All assemblies complete! Results:"
echo "======================================"
cat assembly_summary.txt
