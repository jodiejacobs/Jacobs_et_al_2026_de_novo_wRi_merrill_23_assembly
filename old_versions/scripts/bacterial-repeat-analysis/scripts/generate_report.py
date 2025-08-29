#!/usr/bin/env python3
"""
Generate a comprehensive HTML report summarizing repetitive elements
found in bacterial genomes from multiple analysis tools.
"""

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import numpy as np
from jinja2 import Template

def parse_args():
    parser = argparse.ArgumentParser(description='Generate repeat analysis report')
    parser.add_argument('--assembly', required=True, help='Input genome assembly')
    parser.add_argument('--repeatmasker', required=True, help='RepeatMasker output file')
    # Updated arguments for ISEScan output
    parser.add_argument('--isescan_tsv', required=True, help='ISEScan TSV output file')
    parser.add_argument('--isescan_gff', required=True, help='ISEScan GFF output file')
    parser.add_argument('--red', required=True, help='RED repeat file')
    parser.add_argument('--interspersed', required=True, help='Interspersed repeats GFF file')
    parser.add_argument('--output', required=True, help='Output HTML report')
    return parser.parse_args()

def parse_assembly_stats(assembly_file):
    """Get basic statistics about the assembly"""
    stats = {'Total Length': 0, 'Num Contigs': 0, 'GC Content': 0, 'N50': 0}
    
    # Read assembly sequences
    sequences = list(SeqIO.parse(assembly_file, "fasta"))
    
    # Calculate statistics
    stats['Num Contigs'] = len(sequences)
    contig_lengths = [len(seq) for seq in sequences]
    stats['Total Length'] = sum(contig_lengths)
    
    # Calculate GC content
    gc_count = sum(seq.seq.count('G') + seq.seq.count('C') for seq in sequences)
    stats['GC Content'] = round(gc_count / stats['Total Length'] * 100, 2)
    
    # Calculate N50
    sorted_lengths = sorted(contig_lengths, reverse=True)
    cumsum = np.cumsum(sorted_lengths)
    n50_idx = np.where(cumsum >= stats['Total Length'] * 0.5)[0][0]
    stats['N50'] = sorted_lengths[n50_idx]
    
    return stats

def parse_repeatmasker(rm_file):
    """Parse RepeatMasker output file"""
    repeats = []
    
    # Skip header lines
    skip_lines = 3
    
    with open(rm_file, 'r') as f:
        lines = f.readlines()[skip_lines:]
        
        for line in lines:
            # Skip empty lines
            if line.strip() == '':
                continue
                
            # Parse RepeatMasker output line
            fields = line.strip().split()
            if len(fields) >= 15:
                repeat = {
                    'score': fields[0],
                    'div': fields[1],
                    'del': fields[2],
                    'ins': fields[3],
                    'query': fields[4],
                    'q_start': int(fields[5]),
                    'q_end': int(fields[6]),
                    'q_left': fields[7],
                    'strand': fields[8],
                    'repeat': fields[9],
                    'r_class': fields[10],
                    'r_start': int(fields[11]),
                    'r_end': int(fields[12]),
                    'r_left': fields[13],
                    'id': fields[14] if len(fields) > 14 else ''
                }
                repeats.append(repeat)
    
    return repeats

def parse_isescan_tsv(isescan_tsv):
    """Parse ISEScan TSV output file"""
    insertion_sequences = []
    
    with open(isescan_tsv, 'r') as f:
        lines = f.readlines()
        
        # Skip header row if present
        if lines and "seqID" in lines[0]:
            lines = lines[1:]
        
        for line in lines:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) >= 5:  # Need at least seqID, family, isBegin, isEnd
                try:
                    # Based on the actual format: seqID family cluster isBegin isEnd ...
                    is_element = {
                        'seqid': fields[0],         # seqID column
                        'family': fields[1],        # family column
                        'cluster': fields[2],       # cluster column
                        'start': int(fields[3]),    # isBegin column
                        'end': int(fields[4]),      # isEnd column
                        'strand': fields[16] if len(fields) > 16 else '+',  # strand column
                        'score': float(fields[20]) if len(fields) > 20 and fields[20] != '-1' else 0.0,  # Using 'type' as a proxy for score
                        'irlen': int(fields[13]) if len(fields) > 13 and fields[13] != '-' else 0,  # irLen column
                        'irpos': fields[22] if len(fields) > 22 else '-',  # tir column
                        'type': 'insertion_sequence'
                    }
                    insertion_sequences.append(is_element)
                except (ValueError, IndexError) as e:
                    print(f"Warning: Could not parse line in ISEScan TSV: {line.strip()}")
                    print(f"Error: {e}")
                    continue
    
    return insertion_sequences

def parse_isescan_gff(isescan_gff):
    """Parse ISEScan GFF output file for additional features"""
    insertion_sequences = []
    
    with open(isescan_gff, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) >= 8:
                # Only process if this is an insertion sequence feature
                if fields[2] == "insertion_sequence":
                    # Extract attributes
                    attrs = {}
                    if len(fields) > 8:
                        attr_parts = fields[8].split(';')
                        for part in attr_parts:
                            if '=' in part:
                                key, value = part.split('=', 1)
                                attrs[key] = value
                    
                    is_element = {
                        'seqid': fields[0],
                        'family': attrs.get('family', 'unknown'),
                        'start': int(fields[3]),
                        'end': int(fields[4]),
                        'strand': fields[6],
                        'score': float(fields[5]) if fields[5] != '.' else 0.0,
                        'irlen': int(attrs.get('ir_length', 0)),
                        'irpos': attrs.get('ir_positions', '-'),
                        'type': 'GFF'
                    }
                    insertion_sequences.append(is_element)
    
    return insertion_sequences

def parse_red(red_file):
    """Parse RED (Repeat Explorer) output file"""
    repeats = []
    
    with open(red_file, 'r') as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) >= 6:
                repeat = {
                    'repeat_id': fields[0],
                    'seqid': fields[1],
                    'start': int(fields[2]),
                    'end': int(fields[3]),
                    'length': int(fields[4]),
                    'score': float(fields[5])
                }
                repeats.append(repeat)
    
    return repeats

def parse_interspersed(gff_file):
    """Parse interspersed repeats GFF file"""
    repeats = []
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) >= 9:
                attrs = {}
                for attr in fields[8].split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attrs[key] = value
                
                repeat = {
                    'seqid': fields[0],
                    'source': fields[1],
                    'type': fields[2],
                    'start': int(fields[3]),
                    'end': int(fields[4]),
                    'score': float(fields[5]),
                    'strand': fields[6],
                    'phase': fields[7],
                    'id': attrs.get('ID', ''),
                    'pair': attrs.get('Pair', ''),
                    'length': int(attrs.get('Length', 0)),
                    'identity': float(attrs.get('Identity', 0))
                }
                repeats.append(repeat)
    
    return repeats

def summarize_repeats(rm_repeats, is_elements, red_repeats, interspersed_repeats, total_length):
    """Summarize repeat findings across tools"""
    summary = {}
    
    # RepeatMasker summary
    if rm_repeats:
        rm_total_bp = sum(r['q_end'] - r['q_start'] + 1 for r in rm_repeats)
        summary['RepeatMasker'] = {
            'count': len(rm_repeats),
            'total_bp': rm_total_bp,
            'percent_genome': round(rm_total_bp / total_length * 100, 2)
        }
    else:
        summary['RepeatMasker'] = {'count': 0, 'total_bp': 0, 'percent_genome': 0}
    
    # ISEScan summary
    if is_elements:
        is_total_bp = sum(e['end'] - e['start'] + 1 for e in is_elements)
        summary['ISEScan'] = {
            'count': len(is_elements),
            'total_bp': is_total_bp,
            'percent_genome': round(is_total_bp / total_length * 100, 2)
        }
    else:
        summary['ISEScan'] = {'count': 0, 'total_bp': 0, 'percent_genome': 0}
    
    # RED summary
    if red_repeats:
        red_total_bp = sum(r['length'] for r in red_repeats)
        summary['RED'] = {
            'count': len(red_repeats),
            'total_bp': red_total_bp,
            'percent_genome': round(red_total_bp / total_length * 100, 2)
        }
    else:
        summary['RED'] = {'count': 0, 'total_bp': 0, 'percent_genome': 0}
    
    # Interspersed repeats summary
    if interspersed_repeats:
        int_total_bp = sum(r['end'] - r['start'] + 1 for r in interspersed_repeats)
        summary['Interspersed'] = {
            'count': len(interspersed_repeats) // 2,  # Each repeat is counted twice (pair)
            'total_bp': int_total_bp,
            'percent_genome': round(int_total_bp / total_length * 100, 2)
        }
    else:
        summary['Interspersed'] = {'count': 0, 'total_bp': 0, 'percent_genome': 0}
    
    return summary

def generate_plots(summary, rm_repeats, is_elements, outdir):
    """Generate plots for the report"""
    plots = {}
    
    # Create output directory for plots
    os.makedirs(outdir, exist_ok=True)
    
    # Summary bar chart
    tool_names = list(summary.keys())
    percentages = [summary[tool]['percent_genome'] for tool in tool_names]
    
    plt.figure(figsize=(10, 6))
    sns.barplot(x=tool_names, y=percentages)
    plt.title('Percentage of Genome Covered by Repetitive Elements')
    plt.xlabel('Analysis Tool')
    plt.ylabel('Genome Coverage (%)')
    plt.tight_layout()
    plot_path = os.path.join(outdir, 'genome_coverage.png')
    plt.savefig(plot_path)
    plots['genome_coverage'] = os.path.basename(plot_path)
    
    # RepeatMasker repeat class distribution
    if rm_repeats:
        rm_classes = [r['r_class'] for r in rm_repeats]
        class_counts = pd.Series(rm_classes).value_counts()
        
        plt.figure(figsize=(10, 6))
        sns.barplot(x=class_counts.index, y=class_counts.values)
        plt.title('RepeatMasker Repeat Class Distribution')
        plt.xlabel('Repeat Class')
        plt.ylabel('Count')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plot_path = os.path.join(outdir, 'repeatmasker_classes.png')
        plt.savefig(plot_path)
        plots['repeatmasker_classes'] = os.path.basename(plot_path)
    
    # ISEScan family distribution
    if is_elements:
        is_families = [e['family'] for e in is_elements]
        family_counts = pd.Series(is_families).value_counts()
        
        plt.figure(figsize=(10, 6))
        sns.barplot(x=family_counts.index, y=family_counts.values)
        plt.title('ISEScan Insertion Sequence Family Distribution')
        plt.xlabel('IS Family')
        plt.ylabel('Count')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plot_path = os.path.join(outdir, 'isescan_families.png')
        plt.savefig(plot_path)
        plots['isescan_families'] = os.path.basename(plot_path)
    
    return plots

def generate_html_report(stats, summary, plots, output_file):
    """Generate HTML report using Jinja2 template"""
    template_str = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Bacterial Genome Repeat Analysis Report</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 20px; }
            h1, h2, h3 { color: #2c3e50; }
            table { border-collapse: collapse; width: 100%; margin-bottom: 20px; }
            th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
            th { background-color: #f2f2f2; }
            tr:nth-child(even) { background-color: #f9f9f9; }
            .summary { display: flex; flex-wrap: wrap; }
            .summary-box { 
                flex: 1; 
                min-width: 200px; 
                margin: 10px; 
                padding: 15px; 
                border-radius: 5px; 
                background-color: #f8f9fa;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }
            .plot-container { margin: 20px 0; }
            .plot-container img { max-width: 100%; height: auto; }
        </style>
    </head>
    <body>
        <h1>Bacterial Genome Repeat Analysis Report</h1>
        
        <h2>Assembly Statistics</h2>
        <table>
            <tr><th>Metric</th><th>Value</th></tr>
            {% for key, value in stats.items() %}
            <tr><td>{{ key }}</td><td>{{ value }}</td></tr>
            {% endfor %}
        </table>
        
        <h2>Repeat Element Summary</h2>
        <div class="summary">
            {% for tool, data in summary.items() %}
            <div class="summary-box">
                <h3>{{ tool }}</h3>
                <p>Number of elements: <strong>{{ data.count }}</strong></p>
                <p>Total base pairs: <strong>{{ data.total_bp }}</strong></p>
                <p>Genome coverage: <strong>{{ data.percent_genome }}%</strong></p>
            </div>
            {% endfor %}
        </div>
        
        <h2>Visualizations</h2>
        {% for plot_id, plot_file in plots.items() %}
        <div class="plot-container">
            <h3>{{ plot_id.replace('_', ' ').title() }}</h3>
            <img src="{{ plot_file }}" alt="{{ plot_id }}">
        </div>
        {% endfor %}
        
        <h2>Tool Details</h2>
        <h3>RepeatMasker</h3>
        <p>RepeatMasker was used to identify known repetitive elements in the genome.</p>
        
        <h3>ISEScan</h3>
        <p>ISEScan identified insertion sequences (IS) and mobile genetic elements.</p>
        
        <h3>RED (Repeat Explorer)</h3>
        <p>RED was used for de novo identification of repetitive elements.</p>
        
        <h3>Interspersed Repeat Analysis</h3>
        <p>Custom analysis was performed to identify interspersed repeats based on self-alignment.</p>
        
        <p><em>Report generated on: {{ date }}</em></p>
    </body>
    </html>
    """
    
    from datetime import datetime
    
    # Generate report
    template = Template(template_str)
    html_content = template.render(
        stats=stats,
        summary=summary,
        plots=plots,
        date=datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    )
    
    # Write report to file
    with open(output_file, 'w') as f:
        f.write(html_content)

def main():
    args = parse_args()
    
    # Parse assembly statistics
    assembly_stats = parse_assembly_stats(args.assembly)
    
    # Parse tool outputs
    rm_repeats = parse_repeatmasker(args.repeatmasker)
    
    # Parse both ISEScan files and combine results
    is_elements_tsv = parse_isescan_tsv(args.isescan_tsv)
    is_elements_gff = parse_isescan_gff(args.isescan_gff)
    
    # Combine and deduplicate ISEScan elements (might have overlap between TSV and GFF)
    is_element_positions = set()
    is_elements = []
    
    for element in is_elements_tsv + is_elements_gff:
        position_key = (element['seqid'], element['start'], element['end'], element['strand'])
        if position_key not in is_element_positions:
            is_element_positions.add(position_key)
            is_elements.append(element)
    
    red_repeats = parse_red(args.red)
    interspersed_repeats = parse_interspersed(args.interspersed)
    
    # Summarize findings
    summary = summarize_repeats(
        rm_repeats, 
        is_elements, 
        red_repeats, 
        interspersed_repeats,
        assembly_stats['Total Length']
    )
    
    # Generate plots
    plots_dir = os.path.dirname(args.output)
    plots = generate_plots(summary, rm_repeats, is_elements, plots_dir)
    
    # Generate HTML report
    generate_html_report(assembly_stats, summary, plots, args.output)
    
    print(f"Report generated at {args.output}")

if __name__ == "__main__":
    main()