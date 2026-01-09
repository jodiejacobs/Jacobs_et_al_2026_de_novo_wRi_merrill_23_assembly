#!/bin/bash
#SBATCH --job-name=multi_assembly
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=assembly_%j.log
#SBATCH --partition=long
set -e

# Parameters
READS="/private/groups/russelllab/jodie/Jacobs_et_al_2026_de_novo_wRi_merrill_23_assembly/data/basecalled/20251218_Merrill23_wRi.wRiM23.fastq.gz"
GENOME_SIZE="1.3m"
MIN_SIZE=1200000  # 1.2 Mb
MAX_SIZE=1400000  # 1.4 Mb
THREADS=32
WORKDIR="multi_assembly_comparison"

# Create working directory
mkdir -p ${WORKDIR}
cd ${WORKDIR}

# Function to check assembly quality
check_assembly() {
    local fasta=$1
    local name=$2
    
    echo "=== Checking ${name} assembly ===" >> ../assembly_summary.txt
    
    # Count contigs
    n_contigs=$(grep -c ">" ${fasta})
    echo "Number of contigs: ${n_contigs}" >> ../assembly_summary.txt
    
    # Get contig sizes
    seqkit fx2tab -l -n ${fasta} | while read contig size; do
        echo "  ${contig}: ${size} bp" >> ../assembly_summary.txt
        
        # Check if it's in the right size range
        if [ ${size} -ge ${MIN_SIZE} ] && [ ${size} -le ${MAX_SIZE} ]; then
            echo "  *** CANDIDATE CHROMOSOME: ${contig} (${size} bp) ***" >> ../assembly_summary.txt
        fi
    done
    
    # Calculate N50
    n50=$(seqkit stats -T ${fasta} | tail -n 1 | cut -f14)
    total=$(seqkit stats -T ${fasta} | tail -n 1 | cut -f5)
    echo "Total size: ${total}" >> ../assembly_summary.txt
    echo "N50: ${n50}" >> ../assembly_summary.txt
    echo "" >> ../assembly_summary.txt
}

# Initialize summary file
echo "Assembly Comparison for wRi Merrill 23" > assembly_summary.txt
echo "Date: $(date)" >> assembly_summary.txt
echo "Reads: ${READS}" >> assembly_summary.txt
echo "Expected chromosome size: ${MIN_SIZE}-${MAX_SIZE} bp" >> assembly_summary.txt
echo "======================================" >> assembly_summary.txt
echo "" >> assembly_summary.txt

##################################################
# 1. Flye (default parameters)
##################################################
echo "Running Flye (default)..."
mkdir -p flye_default
mamba run -n assembly flye \
    --nano-hq ${READS} \
    --genome-size ${GENOME_SIZE} \
    --threads ${THREADS} \
    --out-dir flye_default

check_assembly flye_default/assembly.fasta "Flye_default"

##################################################
# 2. Flye (conservative - higher min-overlap)
##################################################
echo "Running Flye (conservative)..."
mkdir -p flye_conservative
mamba run -n assembly flye \
    --nano-hq ${READS} \
    --genome-size ${GENOME_SIZE} \
    --threads ${THREADS} \
    --min-overlap 5000 \
    --out-dir flye_conservative

check_assembly flye_conservative/assembly.fasta "Flye_conservative"

##################################################
# 3. Canu
##################################################
echo "Running Canu..."
mkdir -p canu_assembly
mamba run -n assembly canu \
    -p wRi -d canu_assembly \
    genomeSize=${GENOME_SIZE} \
    -nanopore ${READS} \
    maxThreads=${THREADS} \
    useGrid=false

check_assembly canu_assembly/wRi.contigs.fasta "Canu"

##################################################
# 4. Raven
##################################################
echo "Running Raven..."
mkdir -p raven_assembly
mamba run -n assembly raven \
    --threads ${THREADS} \
    ${READS} > raven_assembly/assembly.fasta

check_assembly raven_assembly/assembly.fasta "Raven"

##################################################
# 5. Miniasm + Minipolish
##################################################
echo "Running Miniasm + Minipolish..."
mkdir -p miniasm_assembly

# Minimap2 overlap
minimap2 -x ava-ont -t ${THREADS} ${READS} ${READS} | \
    gzip -1 > miniasm_assembly/overlaps.paf.gz

# Miniasm
miniasm -f ${READS} miniasm_assembly/overlaps.paf.gz > miniasm_assembly/assembly.gfa

# Convert GFA to FASTA
awk '/^S/{print ">"$2; print $3}' miniasm_assembly/assembly.gfa > miniasm_assembly/assembly_unpolished.fasta

# Minipolish
mamba run -n assembly minipolish \
    --threads ${THREADS} \
    ${READS} miniasm_assembly/assembly.gfa > miniasm_assembly/assembly.fasta

check_assembly miniasm_assembly/assembly.fasta "Miniasm+Minipolish"

##################################################
# 6. Hifiasm (ONT mode)
##################################################
echo "Running Hifiasm (ONT mode)..."
mkdir -p hifiasm_assembly
cd hifiasm_assembly

mamba run -n assembly hifiasm \
    -o wRi \
    -t ${THREADS} \
    --ul ../${READS}

# Convert GFA to FASTA
awk '/^S/{print ">"$2; print $3}' wRi.bp.p_ctg.gfa > assembly.fasta

cd ..
check_assembly hifiasm_assembly/assembly.fasta "Hifiasm"

##################################################
# Check for circular contigs using minimap2 self-alignment
##################################################
echo "" >> assembly_summary.txt
echo "=== CHECKING FOR CIRCULAR CONTIGS ===" >> assembly_summary.txt
echo "" >> assembly_summary.txt

for asm_dir in flye_default flye_conservative canu_assembly raven_assembly miniasm_assembly hifiasm_assembly; do
    if [ -d ${asm_dir} ]; then
        # Find the assembly file
        if [ -f ${asm_dir}/assembly.fasta ]; then
            asm_file=${asm_dir}/assembly.fasta
        elif [ -f ${asm_dir}/wRi.contigs.fasta ]; then
            asm_file=${asm_dir}/wRi.contigs.fasta
        else
            continue
        fi
        
        echo "Checking circularity: ${asm_dir}" >> ../assembly_summary.txt
        
        # For each contig in size range, check if ends overlap
        seqkit fx2tab -l -n ${asm_file} | while read contig size; do
            if [ ${size} -ge ${MIN_SIZE} ] && [ ${size} -le ${MAX_SIZE} ]; then
                echo "  Checking ${contig} (${size} bp)..." >> ../assembly_summary.txt
                
                # Extract first and last 5kb
                seqkit subseq -r 1:5000 ${asm_file} -j ${THREADS} --id-regexp "^${contig}$" > ${asm_dir}/${contig}_start.fa 2>/dev/null || continue
                seqkit subseq -r -5000:-1 ${asm_file} -j ${THREADS} --id-regexp "^${contig}$" > ${asm_dir}/${contig}_end.fa 2>/dev/null || continue
                
                # Check for overlap
                overlap=$(minimap2 -x asm5 ${asm_dir}/${contig}_start.fa ${asm_dir}/${contig}_end.fa 2>/dev/null | wc -l)
                
                if [ ${overlap} -gt 0 ]; then
                    echo "    *** CIRCULAR: Detected ${overlap} overlapping alignment(s) ***" >> ../assembly_summary.txt
                else
                    echo "    No overlap detected (may still be circular)" >> ../assembly_summary.txt
                fi
            fi
        done
        echo "" >> ../assembly_summary.txt
    fi
done

##################################################
# Create comparison table
##################################################
echo "" >> assembly_summary.txt
echo "=== SUMMARY TABLE ===" >> assembly_summary.txt
echo "Assembler | N_Contigs | Total_Size | N50 | Candidate_Chromosomes" >> assembly_summary.txt
echo "----------|-----------|------------|-----|----------------------" >> assembly_summary.txt

for asm_dir in flye_default flye_conservative canu_assembly raven_assembly miniasm_assembly hifiasm_assembly; do
    if [ -d ${asm_dir} ]; then
        if [ -f ${asm_dir}/assembly.fasta ]; then
            asm_file=${asm_dir}/assembly.fasta
        elif [ -f ${asm_dir}/wRi.contigs.fasta ]; then
            asm_file=${asm_dir}/wRi.contigs.fasta
        else
            continue
        fi
        
        n_contigs=$(grep -c ">" ${asm_file})
        stats=$(seqkit stats -T ${asm_file} | tail -n 1)
        total=$(echo ${stats} | cut -f5)
        n50=$(echo ${stats} | cut -f14)
        
        # Count candidates
        candidates=$(seqkit fx2tab -l -n ${asm_file} | awk -v min=${MIN_SIZE} -v max=${MAX_SIZE} '$2 >= min && $2 <= max' | wc -l)
        
        echo "${asm_dir} | ${n_contigs} | ${total} | ${n50} | ${candidates}" >> ../assembly_summary.txt
    fi
done

echo ""
echo "Assembly complete! Check assembly_summary.txt for results."
cat assembly_summary.txt
