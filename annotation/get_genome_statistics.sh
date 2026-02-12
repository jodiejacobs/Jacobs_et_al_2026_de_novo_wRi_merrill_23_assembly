#!/bin/bash

FASTA="/Users/jodiejacobs/Library/CloudStorage/GoogleDrive-jomojaco@ucsc.edu/My Drive/wRi_Merrill_23_paper/new_assembly/20251218_Merrill23_wRi_M23.assembly.fasta"
GFF="/Users/jodiejacobs/Library/CloudStorage/GoogleDrive-jomojaco@ucsc.edu/My Drive/wRi_Merrill_23_paper/wRi_Merrill_23_prokka.gff"

echo "=== Genome Statistics ==="
echo ""

# Length
LENGTH=$(grep -v "^>" "$FASTA" | tr -d '\n' | wc -c | tr -d ' ')
echo "Length (bp): $(printf "%'d" $LENGTH)"

# GC content
SEQ=$(grep -v "^>" "$FASTA" | tr -d '\n' | tr '[:lower:]' '[:upper:]')
GC=$(echo "$SEQ" | grep -o "[GC]" | wc -l | tr -d ' ')
TOTAL=$(echo "$SEQ" | wc -c | tr -d ' ')
GC_PERCENT=$(echo "scale=2; $GC * 100 / $TOTAL" | bc)
echo "GC Content: ${GC_PERCENT}%"

# From GFF - count features by type
CDS=$(grep -v '^#' "$GFF" | awk '$3=="CDS"' | wc -l | tr -d ' ')
TRNA=$(grep -v '^#' "$GFF" | awk '$3=="tRNA"' | wc -l | tr -d ' ')
RRNA=$(grep -v '^#' "$GFF" | awk '$3=="rRNA"' | wc -l | tr -d ' ')
NCRNA=$(grep -v '^#' "$GFF" | awk '$3=="ncRNA"' | wc -l | tr -d ' ')

# Pseudogenes (look for /pseudo in attributes or specific product names)
PSEUDO=$(grep -v '^#' "$GFF" | grep -i "pseudogene\|frameshifted\|truncat" | wc -l | tr -d ' ')

# Total genes = all protein-coding + RNA genes
# For prokaryotes, this is typically CDS + tRNA + rRNA + ncRNA
TOTAL_GENES=$((CDS + TRNA + RRNA + NCRNA))
RNA_TOTAL=$((RRNA + TRNA + NCRNA))

echo "Genes (total): $TOTAL_GENES"
echo "CDSs (total): $CDS"
echo "Genes (RNA): $RNA_TOTAL"

# rRNA breakdown - count by searching product annotation
S5=$(grep -v '^#' "$GFF" | awk '$3=="rRNA"' | grep -c "5S")
S16=$(grep -v '^#' "$GFF" | awk '$3=="rRNA"' | grep -c "16S")
S23=$(grep -v '^#' "$GFF" | awk '$3=="rRNA"' | grep -c "23S")
echo "rRNAs: $S5, $S16, $S23 (5S, 16S, 23S)"

echo "tRNAs: $TRNA"
echo "ncRNAs: $NCRNA"
echo "Pseudogenes (total): $PSEUDO"
