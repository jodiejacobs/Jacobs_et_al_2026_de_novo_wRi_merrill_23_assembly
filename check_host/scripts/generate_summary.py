# scripts/generate_summary.py

import pandas as pd
import sys
from pathlib import Path

# Parse inputs
sample = snakemake.params.sample
idxstats_file = snakemake.input.idxstats
organisms = snakemake.params.organisms
ref_dict = snakemake.params.ref_dict
avg_read_length = snakemake.params.avg_read_length

# Read idxstats
data = []
with open(idxstats_file) as f:
    for line in f:
        if line.startswith('*'):
            continue
        fields = line.strip().split('\t')
        ref = fields[0]
        length = int(fields[1])
        mapped = int(fields[2])
        unmapped = int(fields[3])
        
        # Extract organism prefix
        if '_' in ref:
            organism = ref.split('_')[0]
        else:
            organism = 'unknown'
        
        data.append({
            'organism': organism,
            'contig': ref,
            'length': length,
            'mapped_reads': mapped,
            'unmapped_reads': unmapped
        })

df = pd.DataFrame(data)

# Summarize by organism
summary = df.groupby('organism').agg({
    'mapped_reads': 'sum',
    'unmapped_reads': 'sum',
    'length': 'sum'
}).reset_index()

# Add descriptions from config
summary['description'] = summary['organism'].map(
    lambda x: ref_dict.get(x, {}).get('description', 'Unknown')
)

# Calculate statistics
total_mapped = summary['mapped_reads'].sum()
summary['percentage'] = 100 * summary['mapped_reads'] / total_mapped
summary['coverage'] = (summary['mapped_reads'] * avg_read_length) / summary['length']

# Sort by mapped reads
summary = summary.sort_values('mapped_reads', ascending=False)

# Write summary report
with open(snakemake.output.summary, 'w') as f:
    f.write(f"Identity Check Summary for {sample}\n")
    f.write("=" * 80 + "\n\n")
    
    f.write("Mapping Statistics:\n")
    f.write("-" * 80 + "\n")
    for _, row in summary.iterrows():
        f.write(f"{row['organism']:10s} ({row['description']:30s}): "
                f"{row['mapped_reads']:>10,} reads ({row['percentage']:>5.2f}%) "
                f"| {row['coverage']:>6.1f}x coverage\n")
    
    f.write("\n" + "=" * 80 + "\n")
    f.write(f"Total mapped reads: {total_mapped:,}\n")
    
    # Interpretation
    f.write("\nInterpretation:\n")
    f.write("-" * 80 + "\n")
    
    # Find most abundant host and Wolbachia
    host_orgs = summary[summary['organism'].str.contains('d[a-z]+', regex=True)]
    wolb_orgs = summary[summary['organism'].str.startswith('w')]
    
    if not host_orgs.empty:
        top_host = host_orgs.iloc[0]
        f.write(f"Primary host: {top_host['organism']} ({top_host['description']}) "
                f"- {top_host['percentage']:.1f}% of reads\n")
        
        # Check for other hosts
        if len(host_orgs) > 1:
            f.write(f"  Alternative hosts detected:\n")
            for _, row in host_orgs.iloc[1:].iterrows():
                if row['percentage'] > 0.1:  # Only report if >0.1%
                    f.write(f"    - {row['organism']} ({row['description']}): "
                            f"{row['percentage']:.2f}%\n")
    
    if not wolb_orgs.empty:
        top_wolb = wolb_orgs.iloc[0]
        f.write(f"Primary Wolbachia: {top_wolb['organism']} ({top_wolb['description']}) "
                f"- {top_wolb['percentage']:.1f}% of reads\n")
        
        # Calculate titer (Wolbachia/Host ratio)
        if not host_orgs.empty:
            titer = top_wolb['mapped_reads'] / top_host['mapped_reads']
            f.write(f"Estimated titer: {titer:.4f} "
                    f"({titer*100:.2f}% Wolbachia/Host)\n")
        
        # Check for other Wolbachia strains
        if len(wolb_orgs) > 1:
            f.write(f"  Alternative Wolbachia strains detected:\n")
            for _, row in wolb_orgs.iloc[1:].iterrows():
                if row['percentage'] > 0.1:  # Only report if >0.1%
                    f.write(f"    - {row['organism']} ({row['description']}): "
                            f"{row['percentage']:.2f}%\n")
    
    # Conclusion
    f.write("\n" + "=" * 80 + "\n")
    f.write("CONCLUSION:\n")
    if not host_orgs.empty and not wolb_orgs.empty:
        f.write(f"Sample appears to be {top_host['organism']} infected with "
                f"{top_wolb['organism']}\n")
    elif not host_orgs.empty:
        f.write(f"Sample appears to be {top_host['organism']} (no Wolbachia detected)\n")
    else:
        f.write("WARNING: Could not determine host/Wolbachia identity\n")

# Write TSV for plotting
summary.to_csv(snakemake.output.tsv, sep='\t', index=False)