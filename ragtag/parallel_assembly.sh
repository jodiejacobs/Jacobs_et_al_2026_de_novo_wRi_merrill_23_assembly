#!/bin/bash
#SBATCH --job-name=multi_assembly
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-6
#SBATCH --output=assembly_%A_%a.log
#SBATCH --partition=medium
set -e

# Parameters
READS="/private/groups/russelllab/jodie/Jacobs_et_al_2026_de_novo_wRi_merrill_23_assembly/data/basecalled/20251218_Merrill23_wRi.wRiM23.fastq.gz"
GENOME_SIZE="1.3m"
MIN_SIZE=1200000  # 1.2 Mb
MAX_SIZE=1400000  # 1.4 Mb
THREADS=8
WORKDIR="multi_assembly_comparison"

# Create working directory
mkdir -p ${WORKDIR}
cd ${WORKDIR}

# Function to check assembly quality
check_assembly() {
    local fasta=$1
    local name=$2
    local outfile=$3
    
    echo "=== Checking ${name} assembly ===" >> ${outfile}
    
    # Count contigs
    n_contigs=$(grep -c ">" ${fasta})
    echo "Number of contigs: ${n_contigs}" >> ${outfile}
    
    # Get contig sizes
    seqkit fx2tab -l -n ${fasta} | while read contig size; do
        echo "  ${contig}: ${size} bp" >> ${outfile}
        
        # Check if it's in the right size range
        if [ ${size} -ge ${MIN_SIZE} ] && [ ${size} -le ${MAX_SIZE} ]; then
            echo "  *** CANDIDATE CHROMOSOME: ${contig} (${size} bp) ***" >> ${outfile}
        fi
    done
    
    # Calculate N50
    n50=$(seqkit stats -T ${fasta} | tail -n 1 | cut -f14)
    total=$(seqkit stats -T ${fasta} | tail -n 1 | cut -f5)
    echo "Total size: ${total}" >> ${outfile}
    echo "N50: ${n50}" >> ${outfile}
    echo "" >> ${outfile}
}

# Function to check circularity
check_circularity() {
    local asm_file=$1
    local asm_dir=$2
    local outfile=$3
    
    echo "Checking circularity: ${asm_dir}" >> ${outfile}
    
    # For each contig in size range, check if ends overlap
    seqkit fx2tab -l -n ${asm_file} | while read contig size; do
        if [ ${size} -ge ${MIN_SIZE} ] && [ ${size} -le ${MAX_SIZE} ]; then
            echo "  Checking ${contig} (${size} bp)..." >> ${outfile}
            
            # Extract first and last 5kb
            seqkit subseq -r 1:5000 ${asm_file} -j ${THREADS} --id-regexp "^${contig}$" > ${asm_dir}/${contig}_start.fa 2>/dev/null || continue
            seqkit subseq -r -5000:-1 ${asm_file} -j ${THREADS} --id-regexp "^${contig}$" > ${asm_dir}/${contig}_end.fa 2>/dev/null || continue
            
            # Check for overlap
            overlap=$(minimap2 -x asm5 ${asm_dir}/${contig}_start.fa ${asm_dir}/${contig}_end.fa 2>/dev/null | wc -l)
            
            if [ ${overlap} -gt 0 ]; then
                echo "    *** CIRCULAR: Detected ${overlap} overlapping alignment(s) ***" >> ${outfile}
            else
                echo "    No overlap detected (may still be circular)" >> ${outfile}
            fi
        fi
    done
    echo "" >> ${outfile}
}

##################################################
# Run different assemblers based on array task ID
##################################################

case ${SLURM_ARRAY_TASK_ID} in

    1)
        ##################################################
        # Flye (default parameters)
        ##################################################
        echo "Running Flye (default)..."
        mkdir -p flye_default
        mamba run -n assembly flye \
            --nano-hq ${READS} \
            --genome-size ${GENOME_SIZE} \
            --threads ${THREADS} \
            --out-dir flye_default

        check_assembly flye_default/assembly.fasta "Flye_default" flye_default_summary.txt
        check_circularity flye_default/assembly.fasta flye_default flye_default_summary.txt
        ;;

    2)
        ##################################################
        # Flye (conservative - higher min-overlap)
        ##################################################
        echo "Running Flye (conservative)..."
        mkdir -p flye_conservative
        mamba run -n assembly flye \
            --nano-hq ${READS} \
            --genome-size ${GENOME_SIZE} \
            --threads ${THREADS} \
            --min-overlap 5000 \
            --out-dir flye_conservative

        check_assembly flye_conservative/assembly.fasta "Flye_conservative" flye_conservative_summary.txt
        check_circularity flye_conservative/assembly.fasta flye_conservative flye_conservative_summary.txt
        ;;

    3)
        ##################################################
        # Canu
        ##################################################
        echo "Running Canu..."
        mkdir -p canu_assembly
        mamba run -n assembly canu \
            -p wRi -d canu_assembly \
            genomeSize=${GENOME_SIZE} \
            -nanopore ${READS} \
            maxThreads=${THREADS} \
            useGrid=false

        check_assembly canu_assembly/wRi.contigs.fasta "Canu" canu_summary.txt
        check_circularity canu_assembly/wRi.contigs.fasta canu_assembly canu_summary.txt
        ;;

    4)
        ##################################################
        # Raven
        ##################################################
        echo "Running Raven..."
        mkdir -p raven_assembly
        mamba run -n assembly raven \
            --threads ${THREADS} \
            ${READS} > raven_assembly/assembly.fasta

        check_assembly raven_assembly/assembly.fasta "Raven" raven_summary.txt
        check_circularity raven_assembly/assembly.fasta raven_assembly raven_summary.txt
        ;;

    5)
        ##################################################
        # Miniasm + Minipolish
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

        check_assembly miniasm_assembly/assembly.fasta "Miniasm+Minipolish" miniasm_summary.txt
        check_circularity miniasm_assembly/assembly.fasta miniasm_assembly miniasm_summary.txt
        ;;

    6)
        ##################################################
        # Hifiasm (ONT mode)
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
        check_assembly hifiasm_assembly/assembly.fasta "Hifiasm" hifiasm_summary.txt
        check_circularity hifiasm_assembly/assembly.fasta hifiasm_assembly hifiasm_summary.txt
        ;;

esac

echo "Job ${SLURM_ARRAY_TASK_ID} complete!"
