#!/usr/bin/env python3
"""
Generate genome plots visualizing repetitive elements from various analysis tools
(RepeatMasker, ISEScan, RED, etc.) for all contigs in the assembly.
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.colors import LinearSegmentedColormap
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
import matplotlib.patches as mpatches
from matplotlib.path import Path
import math

def parse_args():
    parser = argparse.ArgumentParser(description='Generate genome visualization with repeats for all contigs')
    parser.add_argument('--assembly', required=True, help='Input genome assembly')
    parser.add_argument('--repeatmasker', required=True, help='RepeatMasker output file')
    # Updated arguments for ISEScan output
    parser.add_argument('--isescan_tsv', required=True, help='ISEScan TSV output file')
    parser.add_argument('--isescan_gff', required=True, help='ISEScan GFF output file')
    parser.add_argument('--red', required=True, help='RED repeat file')
    parser.add_argument('--interspersed', required=True, help='Interspersed repeats GFF file')
    parser.add_argument('--output', required=True, help='Output plot file basename')
    parser.add_argument('--density_output', required=True, help='Output density plot file basename')
    parser.add_argument('--min_contig_length', type=int, default=1000, help='Minimum contig length to plot (default: 1kb)')
    parser.add_argument('--max_contigs', type=int, default=10, help='Maximum number of contigs to plot (default: 20)')
    parser.add_argument('--circular_contigs', nargs='+', default=[], help='List of contig IDs that are known to be circular')
    return parser.parse_args()

def plot_circular_genome(sequence, features, output_file, is_alternative=False):
    """Generate circular genome plot with repeats using matplotlib"""
    fig, ax = plt.subplots(figsize=(12, 12), subplot_kw={'polar': True})
    
    # Set up the circular plot
    genome_length = len(sequence)
    ax.set_theta_direction(-1)  # clockwise
    ax.set_theta_zero_location('N')  # 0 radians at the top
    
    # Remove radial labels and grid
    ax.set_yticklabels([])
    ax.grid(False)
    
    # Draw genome backbone
    ax.plot([0, 2*np.pi], [1, 1], color='black', linewidth=1)
    
    # Get colors for feature types
    colors = get_feature_colors()
    
    # Filter features for this specific contig
    contig_features = [f for f in features if f['contig'] == sequence.id]
    feature_count = len(contig_features)
    
    # Plot features
    for feature in contig_features:
        start_angle = 2 * np.pi * feature['start'] / genome_length
        end_angle = 2 * np.pi * feature['end'] / genome_length
        
        # Determine track radius based on feature type to separate different tools
        radius = 1.0
        track_offset = {
            'RepeatMasker': 0.1,
            'ISEScan': 0.2,
            'RED': 0.3,
            'Interspersed': 0.4
        }
        radius_offset = track_offset.get(feature['type'], 0.1)
        
        # Draw feature arc
        arc = patches.Wedge(
            center=(0, 0),
            r=radius + radius_offset,
            theta1=np.degrees(start_angle),
            theta2=np.degrees(end_angle),
            width=0.05,
            color=colors.get(feature['type'], '#CCCCCC')
        )
        ax.add_patch(arc)
    
    # Add ticks for genome position
    tick_positions = np.linspace(0, 2*np.pi, 13)[:-1]  # 12 positions (0, 30, 60, ... degrees)
    tick_labels = [f"{int(i * genome_length / (2*np.pi) / 1000)}k" for i in tick_positions]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    
    # Add legend
    legend_elements = []
    for feature_type, color in colors.items():
        patch = mpatches.Patch(color=color, label=feature_type)
        legend_elements.append(patch)
    
    plt.legend(handles=legend_elements, loc='center', bbox_to_anchor=(0.5, 0.5), ncol=1)
    
    # Set title with appropriate indication
    is_circular = determine_if_circular(sequence)
    if is_circular:
        title_prefix = "Circular Genome Map"
    else:
        title_prefix = "Circular Visualization of Linear Genome"
    
    if is_alternative:
        title_suffix = "\n(alternative visualization)"
    else:
        title_suffix = ""
    
    plt.title(f"{title_prefix} with Repetitive Elements\n{sequence.id}, {len(sequence):,} bp\n{feature_count} features{title_suffix}", 
              fontsize=14, y=1.05)
    
    # Save figure
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    message = "alternative" if is_alternative else "main"
    print(f"Circular genome plot saved to {output_file} ({message})")
    return feature_count

def get_assembly_contigs(assembly_file, min_length=10000, max_contigs=10):
    """Get all contigs from the assembly that meet the minimum length criteria"""
    sequences = list(SeqIO.parse(assembly_file, "fasta"))
    
    # Filter by minimum length
    filtered_sequences = [seq for seq in sequences if len(seq) >= min_length]
    
    # Sort by length (largest first)
    filtered_sequences.sort(key=lambda seq: -len(seq))
    
    # Limit to max_contigs
    if max_contigs > 0 and len(filtered_sequences) > max_contigs:
        filtered_sequences = filtered_sequences[:max_contigs]
        print(f"Limiting visualization to the {max_contigs} largest contigs")
    
    print(f"Visualizing {len(filtered_sequences)} contigs from the assembly")
    return filtered_sequences

def determine_if_circular(sequence, circular_contigs=None):
    """Determine if a contig is circular based on its description, ID, or explicit list"""
    # First check explicit list of circular contigs if provided
    if circular_contigs and sequence.id in circular_contigs:
        return True
    
    # Then look for keywords in the description or ID that suggest circularity
    circular_keywords = ['circular', 'complete', 'plasmid', 'chromosome', 'circle']
    
    # Check the description (if available)
    if hasattr(sequence, 'description'):
        desc_lower = sequence.description.lower()
        for keyword in circular_keywords:
            if keyword in desc_lower:
                return True
    
    # Check the ID
    id_lower = sequence.id.lower()
    for keyword in circular_keywords:
        if keyword in id_lower:
            return True
    
    # Default to linear if no evidence of circularity
    return False

def parse_repeatmasker(rm_file):
    """Parse RepeatMasker output file"""
    features = []
    
    # Skip header lines
    skip_lines = 3
    
    try:
        with open(rm_file, 'r') as f:
            lines = f.readlines()[skip_lines:]
            
            for line in lines:
                # Skip empty lines
                if line.strip() == '':
                    continue
                    
                # Parse RepeatMasker output line
                fields = line.strip().split()
                if len(fields) >= 15:
                    # Extract contig ID from the query name
                    contig_id = fields[4]
                    
                    feature = {
                        'contig': contig_id,
                        'start': int(fields[5]),
                        'end': int(fields[6]),
                        'strand': 1 if fields[8] == '+' else -1,
                        'type': 'RepeatMasker',
                        'subtype': fields[10],  # Repeat class
                        'label': f"RM:{fields[9]}"  # Repeat name
                    }
                    features.append(feature)
    except Exception as e:
        print(f"Warning: Error processing RepeatMasker file: {e}")
    
    return features

def parse_isescan_tsv(isescan_tsv):
    """Parse ISEScan TSV output file"""
    features = []
    
    try:
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
                        # Extract contig ID
                        contig_id = fields[0]
                        
                        # Based on the actual format: seqID family cluster isBegin isEnd ...
                        feature = {
                            'contig': contig_id,
                            'start': int(fields[3]),    # isBegin column
                            'end': int(fields[4]),      # isEnd column
                            'strand': 1 if len(fields) > 16 and fields[16] == '+' else -1 if len(fields) > 16 and fields[16] == '-' else 1,  # strand column
                            'type': 'ISEScan',
                            'subtype': fields[1],       # family column
                            'label': f"IS:{fields[1]}"
                        }
                        features.append(feature)
                    except (ValueError, IndexError) as e:
                        print(f"Warning: Could not parse line in ISEScan TSV: {line.strip()}")
                        print(f"Error: {e}")
                        continue
    except Exception as e:
        print(f"Warning: Error processing ISEScan TSV file: {e}")
    
    return features

def parse_isescan_gff(isescan_gff):
    """Parse ISEScan GFF output file for additional features"""
    features = []
    
    try:
        with open(isescan_gff, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) >= 8:
                    # Only process if this is an insertion sequence feature
                    if fields[2].lower() == "insertion_sequence":
                        # Extract contig ID from first column
                        contig_id = fields[0]
                        
                        # Extract attributes
                        attrs = {}
                        if len(fields) > 8:
                            attr_parts = fields[8].split(';')
                            for part in attr_parts:
                                if '=' in part:
                                    key, value = part.split('=', 1)
                                    attrs[key] = value
                        
                        feature = {
                            'contig': contig_id,
                            'start': int(fields[3]),
                            'end': int(fields[4]),
                            'strand': 1 if fields[6] == '+' else -1,
                            'type': 'ISEScan',
                            'subtype': attrs.get('family', 'unknown'),
                            'label': f"IS:{attrs.get('family', 'unknown')}"
                        }
                        features.append(feature)
    except Exception as e:
        print(f"Warning: Error processing ISEScan GFF file: {e}")
    
    return features

def parse_red(red_file):
    """Parse RED (Repeat Explorer) output file"""
    features = []
    
    try:
        with open(red_file, 'r') as f:
            for line in f:
                fields = line.strip().split()
                if len(fields) >= 6:
                    # Extract contig ID from the sequence ID field
                    contig_id = fields[1]
                    
                    feature = {
                        'contig': contig_id,
                        'start': int(fields[2]),
                        'end': int(fields[3]),
                        'strand': 1,  # RED doesn't provide strand info
                        'type': 'RED',
                        'subtype': 'de_novo',
                        'label': f"RED:{fields[0]}"
                    }
                    features.append(feature)
    except Exception as e:
        print(f"Warning: Error processing RED file: {e}")
    
    return features

def parse_interspersed(gff_file):
    """Parse interspersed repeats GFF file"""
    features = []
    
    try:
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) >= 9:
                    # Extract contig ID from first column
                    contig_id = fields[0]
                    
                    attrs = {}
                    for attr in fields[8].split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            attrs[key] = value
                    
                    feature = {
                        'contig': contig_id,
                        'start': int(fields[3]),
                        'end': int(fields[4]),
                        'strand': 1 if fields[6] == '+' else -1,
                        'type': 'Interspersed',
                        'subtype': 'repeat',
                        'label': f"INT:{attrs.get('ID', '')}"
                    }
                    features.append(feature)
    except Exception as e:
        print(f"Warning: Error processing interspersed repeats file: {e}")
    
    return features

def get_feature_colors():
    """Define colors for different feature types"""
    colors = {
        'RepeatMasker': '#FF5733',  # Orange-red
        'ISEScan': '#33A8FF',       # Blue
        'RED': '#33FF57',           # Green
        'Interspersed': '#A833FF'   # Purple
    }
    return colors

def plot_genome(sequence, features, output_file):
    """Generate appropriate genome plot based on whether the contig is circular or linear"""
    # Determine if the sequence is circular
    is_circular = determine_if_circular(sequence)
    
    if is_circular:
        # Use circular plot for circular contigs
        return plot_circular_genome(sequence, features, output_file)
    else:
        # Use linear plot for linear contigs, but still create both plot types
        feature_count = plot_linear_as_main(sequence, features, output_file)
        # Create circular view as alternative visualization (with "_circular" suffix)
        circular_output = output_file.replace('.pdf', '_circular.pdf')
        plot_circular_genome(sequence, features, circular_output, is_alternative=True)
        return feature_count

def plot_linear_as_main(sequence, features, output_file):
    """Create a linear genome plot as the main visualization"""
    fig, ax = plt.subplots(figsize=(15, 8))
    
    # Get genome length and colors
    genome_length = len(sequence)
    colors = get_feature_colors()
    
    # Filter features for this specific contig
    contig_features = [f for f in features if f['contig'] == sequence.id]
    feature_count = len(contig_features)
    
    # Set up tracks for different feature types
    tracks = {
        'RepeatMasker': 0,
        'ISEScan': 1,
        'RED': 2,
        'Interspersed': 3
    }
    num_tracks = len(tracks)
    
    # Draw genome backbone
    ax.axhline(y=num_tracks + 0.5, color='black', linestyle='-', linewidth=1)
    
    # Plot features
    for feature in contig_features:
        track = tracks.get(feature['type'], 0)
        height = 0.7  # Height of the feature rectangle
        
        # Draw feature as rectangle
        rect = patches.Rectangle(
            (feature['start'], track),
            width=feature['end'] - feature['start'],
            height=height,
            color=colors.get(feature['type'], '#CCCCCC')
        )
        ax.add_patch(rect)
    
    # Set axis limits
    ax.set_xlim(0, genome_length)
    ax.set_ylim(-0.5, num_tracks + 1)
    
    # Add track labels
    track_labels = list(tracks.keys())
    plt.yticks(
        [tracks[label] + 0.35 for label in track_labels],
        track_labels
    )
    
    # Format x-axis with kilobase ticks
    kilobase_ticks = np.linspace(0, genome_length, 10)
    kilobase_labels = [f"{int(pos/1000)}k" for pos in kilobase_ticks]
    plt.xticks(kilobase_ticks, kilobase_labels)
    
    # Add title and labels
    plt.title(f"Linear Genome Map with Repetitive Elements\n{sequence.id}, {len(sequence):,} bp\n{feature_count} features", fontsize=14)
    plt.xlabel('Genome Position (bp)')
    
    # Add legend
    legend_elements = []
    for feature_type, color in colors.items():
        patch = mpatches.Patch(color=color, label=feature_type)
        legend_elements.append(patch)
    
    ax.legend(handles=legend_elements, loc='upper right')
    
    # Save figure
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print(f"Linear genome plot saved to {output_file}")
    return feature_count

def create_linear_plot(sequence, features, output_file):
    """Create a linear genome plot for more detailed view using matplotlib"""
    fig, ax = plt.subplots(figsize=(15, 8))
    
    # Get genome length and colors
    genome_length = len(sequence)
    colors = get_feature_colors()
    
    # Filter features for this specific contig
    contig_features = [f for f in features if f['contig'] == sequence.id]
    feature_count = len(contig_features)
    
    # Set up tracks for different feature types
    tracks = {
        'RepeatMasker': 0,
        'ISEScan': 1,
        'RED': 2,
        'Interspersed': 3
    }
    num_tracks = len(tracks)
    
    # Draw genome backbone
    ax.axhline(y=num_tracks + 0.5, color='black', linestyle='-', linewidth=1)
    
    # Plot features
    for feature in contig_features:
        track = tracks.get(feature['type'], 0)
        height = 0.7  # Height of the feature rectangle
        
        # Draw feature as rectangle
        rect = patches.Rectangle(
            (feature['start'], track),
            width=feature['end'] - feature['start'],
            height=height,
            color=colors.get(feature['type'], '#CCCCCC')
        )
        ax.add_patch(rect)
    
    # Set axis limits
    ax.set_xlim(0, genome_length)
    ax.set_ylim(-0.5, num_tracks + 1)
    
    # Add track labels
    track_labels = list(tracks.keys())
    plt.yticks(
        [tracks[label] + 0.35 for label in track_labels],
        track_labels
    )
    
    # Format x-axis with kilobase ticks
    kilobase_ticks = np.linspace(0, genome_length, 10)
    kilobase_labels = [f"{int(pos/1000)}k" for pos in kilobase_ticks]
    plt.xticks(kilobase_ticks, kilobase_labels)
    
    # Add title and labels
    plt.title(f"Linear Genome Map with Repetitive Elements\n{sequence.id}, {len(sequence):,} bp\n{feature_count} features", fontsize=14)
    plt.xlabel('Genome Position (bp)')
    
    # Save figure
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print(f"Linear genome plot saved to {output_file}")

def create_repeat_density_plot(sequence, features, output_file):
    """Create a repeat density plot to show distribution of repeats"""
    # Set up the figure
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Filter features for this specific contig
    contig_features = [f for f in features if f['contig'] == sequence.id]
    
    # Prepare data for density plot
    seq_length = len(sequence)
    window_size = max(1000, seq_length // 1000)  # Adaptive window size
    windows = np.arange(0, seq_length, window_size)
    
    # Count features in each window
    feature_counts = {}
    colors = get_feature_colors()
    
    for feature_type in colors.keys():
        counts = [0] * len(windows)
        type_features = [f for f in contig_features if f['type'] == feature_type]
        
        for feature in type_features:
            for i, start in enumerate(windows):
                end = start + window_size
                # Check for overlap
                if feature['start'] < end and feature['end'] > start:
                    counts[i] += 1
        
        feature_counts[feature_type] = counts
    
    # Plot the density
    bottom = np.zeros(len(windows))
    for feature_type, counts in feature_counts.items():
        ax.bar(windows, counts, bottom=bottom, width=window_size, 
               label=feature_type, color=colors[feature_type], alpha=0.7)
        bottom += np.array(counts)
    
    # Customize plot
    ax.set_xlabel('Genome Position (bp)')
    ax.set_ylabel('Number of Repeat Features')
    ax.set_title(f'Repeat Density Across the Genome\n{sequence.id}, {len(sequence):,} bp')
    ax.legend()
    
    # Set x-axis limits
    ax.set_xlim(0, seq_length)
    
    # Add genome length ticks
    x_ticks = np.linspace(0, seq_length, 10)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f'{int(x/1000)}k' for x in x_ticks], rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close(fig)
    
    print(f"Density plot saved to {output_file}")

def create_repeat_summary_plot(contigs, all_features, output_file):
    """Create a summary plot showing repeat content across all contigs"""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Get colors for feature types
    colors = get_feature_colors()
    feature_types = list(colors.keys())
    
    # Prepare data for stacked bar chart
    contig_names = [seq.id for seq in contigs]
    contig_lengths = [len(seq) for seq in contigs]
    
    # Create a dictionary to store counts by contig and feature type
    data = {contig_id: {feature_type: 0 for feature_type in feature_types} for contig_id in contig_names}
    
    # Count features by type for each contig
    for feature in all_features:
        contig = feature['contig']
        feature_type = feature['type']
        if contig in data and feature_type in data[contig]:
            # Count base pairs covered by the feature
            feature_length = feature['end'] - feature['start']
            data[contig][feature_type] += feature_length
    
    # Convert counts to percentages of contig length
    for i, contig_id in enumerate(contig_names):
        for feature_type in feature_types:
            data[contig_id][feature_type] = (data[contig_id][feature_type] / contig_lengths[i]) * 100
    
    # Prepare data for plotting
    x = np.arange(len(contig_names))
    width = 0.7
    
    # Create stacked bar chart
    bottom = np.zeros(len(contig_names))
    for feature_type in feature_types:
        values = [data[contig_id][feature_type] for contig_id in contig_names]
        ax.bar(x, values, width, bottom=bottom, label=feature_type, color=colors[feature_type])
        bottom += np.array(values)
    
    # Customize plot
    ax.set_ylabel('Percentage of Contig Covered by Repeats')
    ax.set_title('Repeat Content Across Contigs')
    ax.set_xticks(x)
    
    # Format contig names to include length
    contig_labels = [f"{name}\n({length:,} bp)" for name, length in zip(contig_names, contig_lengths)]
    ax.set_xticklabels(contig_labels, rotation=45, ha='right')
    
    ax.legend()
    
    # Add percentage labels on top of bars
    for i, total in enumerate(bottom):
        ax.text(i, total + 1, f"{total:.1f}%", ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close(fig)
    
    print(f"Summary plot saved to {output_file}")

def main():
    args = parse_args()
    
    # Get all contigs from the assembly
    contigs = get_assembly_contigs(args.assembly, args.min_contig_length, args.max_contigs)
    
    # Get list of explicitly specified circular contigs
    circular_contigs = args.circular_contigs if args.circular_contigs else []
    
    # Parse features from different tools
    print("Parsing RepeatMasker output...")
    rm_features = parse_repeatmasker(args.repeatmasker)
    
    # Parse both ISEScan output files and combine results
    print("Parsing ISEScan TSV output...")
    is_features_tsv = parse_isescan_tsv(args.isescan_tsv)
    print("Parsing ISEScan GFF output...")
    is_features_gff = parse_isescan_gff(args.isescan_gff)
    
    # Combine and deduplicate ISEScan features (might have overlap between TSV and GFF)
    is_feature_positions = set()
    is_features = []
    
    for feature in is_features_tsv + is_features_gff:
        position_key = (feature['contig'], feature['start'], feature['end'], feature['strand'])
        if position_key not in is_feature_positions:
            is_feature_positions.add(position_key)
            is_features.append(feature)
    
    print("Parsing RED output...")
    red_features = parse_red(args.red)
    print("Parsing interspersed repeats...")
    int_features = parse_interspersed(args.interspersed)
    
    # Combine all features
    all_features = rm_features + is_features + red_features + int_features
    print(f"Total features: {len(all_features)}")
    
    # Create output directories
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    os.makedirs(os.path.dirname(args.density_output), exist_ok=True)
    
    # Process each contig separately
    contig_feature_counts = {}
    
    for i, contig in enumerate(contigs):
        print(f"\nProcessing contig {i+1}/{len(contigs)}: {contig.id}")
        
        # Check if contig is circular
        is_circular = determine_if_circular(contig, circular_contigs)
        print(f"  Contig topology: {'circular' if is_circular else 'linear'}")
        
        # Generate appropriate plot for this contig
        if is_circular:
            # For circular contigs, use circular plot as primary
            circular_output = args.output.replace('.pdf', f'_{contig.id}.pdf')
            feature_count = plot_circular_genome(contig, all_features, circular_output)
            contig_feature_counts[contig.id] = feature_count
            
            # Still create linear plot as alternative view
            linear_output = args.output.replace('.pdf', f'_{contig.id}_linear.pdf')
            create_linear_plot(contig, all_features, linear_output)
        else:
            # For linear contigs, use linear plot as primary
            linear_output = args.output.replace('.pdf', f'_{contig.id}.pdf')
            feature_count = plot_linear_as_main(contig, all_features, linear_output)
            contig_feature_counts[contig.id] = feature_count
            
            # Create circular plot as alternative view
            circular_output = args.output.replace('.pdf', f'_{contig.id}_circular.pdf')
            plot_circular_genome(contig, all_features, circular_output, is_alternative=True)
        
        # Generate density plot for this contig
        density_output = args.density_output.replace('.pdf', f'_{contig.id}.pdf')
        create_repeat_density_plot(contig, all_features, density_output)
    
    # Create a summary plot showing all contigs
    summary_output = args.output.replace('.pdf', '_summary.pdf')
    create_repeat_summary_plot(contigs, all_features, summary_output)
    
    # Print a summary of feature counts by contig
    print("\nFeature count summary:")
    for contig in contigs:
        is_circular = determine_if_circular(contig, circular_contigs)
        topology = "circular" if is_circular else "linear"
        print(f"  {contig.id} ({topology}): {contig_feature_counts.get(contig.id, 0)} features")

    print("\nAll plots generated successfully")
    
    print(f"\nAll plots generated successfully")

if __name__ == "__main__":
    main()