#!/usr/bin/env python3
"""
Filter self-BLAST results to identify interspersed repeats in bacterial genomes.
Converts filtered results to GFF format for visualization.
"""

import argparse
import pandas as pd
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description='Filter self-BLAST results to identify repeats')
    parser.add_argument('--input', required=True, help='Input self-BLAST output file')
    parser.add_argument('--output', required=True, help='Output GFF file')
    parser.add_argument('--min_length', type=int, default=50, help='Minimum repeat length')
    parser.add_argument('--min_identity', type=float, default=90, help='Minimum percent identity')
    return parser.parse_args()

def filter_blast_results(blast_file, min_length, min_identity):
    """
    Filter BLAST results for interspersed repeats
    """
    # Read BLAST output
    cols = ['qseqid', 'qstart', 'qend', 'sseqid', 'sstart', 'send', 'pident', 'length']
    df = pd.read_csv(blast_file, sep='\t', names=cols)
    
    # Filter by minimum length and identity
    df = df[(df['length'] >= min_length) & (df['pident'] >= min_identity)]
    
    # Remove self-matches (exact same coordinates)
    df = df[~((df['qseqid'] == df['sseqid']) & 
              (df['qstart'] == df['sstart']) & 
              (df['qend'] == df['send']))]
    
    # Identify interspersed repeats (not tandem)
    # Keep only hits where the distance between query and subject is greater than repeat length
    df['q_distance'] = np.abs(df['qstart'] - df['sstart'])
    df = df[df['q_distance'] > df['length']]
    
    return df

def convert_to_gff(df):
    """
    Convert filtered BLAST results to GFF format
    """
    gff_records = []
    repeat_id = 1
    
    for _, row in df.iterrows():
        # Query feature
        q_record = {
            'seqid': row['qseqid'],
            'source': 'self_blast',
            'type': 'repeat_region',
            'start': int(row['qstart']),
            'end': int(row['qend']),
            'score': row['pident'],
            'strand': '+',
            'phase': '.',
            'attributes': f'ID=repeat_{repeat_id};Pair=repeat_{repeat_id + 1};Length={row["length"]};Identity={row["pident"]}'
        }
        
        # Subject feature
        s_record = {
            'seqid': row['sseqid'],
            'source': 'self_blast',
            'type': 'repeat_region',
            'start': int(row['sstart']),
            'end': int(row['send']),
            'score': row['pident'],
            'strand': '+',
            'phase': '.',
            'attributes': f'ID=repeat_{repeat_id + 1};Pair=repeat_{repeat_id};Length={row["length"]};Identity={row["pident"]}'
        }
        
        gff_records.append(q_record)
        gff_records.append(s_record)
        repeat_id += 2
    
    return gff_records

def write_gff(gff_records, output_file):
    """
    Write GFF records to file
    """
    with open(output_file, 'w') as f:
        # Write GFF3 header
        f.write("##gff-version 3\n")
        
        # Write records
        for record in gff_records:
            line = [
                record['seqid'],
                record['source'],
                record['type'],
                str(record['start']),
                str(record['end']),
                str(record['score']),
                record['strand'],
                record['phase'],
                record['attributes']
            ]
            f.write('\t'.join(line) + '\n')

def main():
    args = parse_args()
    
    # Filter BLAST results
    filtered_df = filter_blast_results(args.input, args.min_length, args.min_identity)
    
    # Convert to GFF
    gff_records = convert_to_gff(filtered_df)
    
    # Write output
    write_gff(gff_records, args.output)
    
    print(f"Found {len(filtered_df)} interspersed repeats. Results written to {args.output}")

if __name__ == "__main__":
    main()
