import os
import sys
import subprocess

def run_command(cmd, description):
    """Runs a shell command and checks for errors."""
    print(f"🔄 Running: {description}...")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"❌ Error: {result.stderr}")
        sys.exit(1)
    return result.stdout

def get_read_positions(bam_file, read_id):
    """Extracts start and end positions of the read from the BAM file."""
    cmd = f"samtools view {bam_file} | grep '{read_id}' | awk '{{print $3, $4, length($10)}}'"
    output = run_command(cmd, "Extracting read positions").strip().split("\n")

    positions = {}
    for line in output:
        fields = line.split()
        if len(fields) >= 3:
            contig, start, read_length = fields[0], int(fields[1]), int(fields[2])
            positions[contig] = (start, start + read_length)

    return positions

def extract_fasta_region(fasta_file, contig, start, end, output_file):
    """Extracts a specific region from a FASTA file."""
    cmd = f"samtools faidx {fasta_file} {contig}:{start}-{end} > {output_file}"
    run_command(cmd, f"Extracting {contig}:{start}-{end}")

def circularize_assembly(assembly_fasta, bridge_fasta, contig_8, contig_18, contig_36, output_fasta):
    # Step 1: Align bridging read to the assembly
    run_command(f"minimap2 -ax map-ont {assembly_fasta} {bridge_fasta} > bridge_read.sam", "Aligning read")
    run_command("samtools view -bS bridge_read.sam | samtools sort -o bridge_read.sorted.bam", "Sorting BAM file")
    run_command("samtools index bridge_read.sorted.bam", "Indexing BAM file")

    # Step 2: Get the alignment positions of the read
    read_id = run_command(f"head -n 1 {bridge_fasta}").strip().replace(">", "")
    positions = get_read_positions("bridge_read.sorted.bam", read_id)

    if contig_8 not in positions or contig_18 not in positions or contig_36 not in positions:
        print("❌ Error: The read does not map to all three contigs.")
        sys.exit(1)

    contig_8_start, contig_8_end = positions[contig_8]
    contig_18_start, contig_18_end = positions[contig_18]
    contig_36_start, contig_36_end = positions[contig_36]

    print(f"✅ Read maps to {contig_8} at {contig_8_start}-{contig_8_end}")
    print(f"✅ Read maps to {contig_18} at {contig_18_start}-{contig_18_end}")
    print(f"✅ Read maps to {contig_36} at {contig_36_start}-{contig_36_end}")

    # Step 3: Extract only the non-overlapping portions of the contigs
    extract_fasta_region(assembly_fasta, contig_8, 1, contig_8_start, "trimmed_contig_8.fasta")
    extract_fasta_region(assembly_fasta, contig_18, contig_18_end, 10000000, "trimmed_contig_18.fasta")
    extract_fasta_region(assembly_fasta, contig_36, contig_36_end, 10000000, "trimmed_contig_36.fasta")

    # Step 4: Extract the non-overlapping part of the bridge read from the FASTA file
    extract_fasta_region(bridge_fasta, read_id, contig_8_end, contig_18_start, "trimmed_bridge_read.fasta")

    # Step 5: Merge all the trimmed sequences correctly
    print("🔗 Merging sequences to form a circularized genome...")
    run_command("cat trimmed_contig_8.fasta trimmed_bridge_read.fasta trimmed_contig_36.fasta trimmed_contig_18.fasta > circularized.fasta",
                "Merging contigs")

    # Step 6: Polish the final assembly
    print("✨ Polishing the circularized genome with Medaka...")
    run_command(f"medaka_consensus -i {bridge_fasta} -d circularized.fasta -o polished_output -m r941_min_hac_g507",
                "Polishing")

    # Step 7: Save the final circularized assembly
    run_command(f"cp polished_output/consensus.fasta {output_fasta}", "Saving final assembly")
    print(f"🎉 Circularized assembly saved as {output_fasta}")

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python circularize_assembly.py assembly.fasta bridge_read.fasta contig_8 contig_18 contig_36 output.fasta")
        sys.exit(1)

    assembly_fasta = sys.argv[1]
    bridge_fasta = sys.argv[2]
    contig_8 = sys.argv[3]
    contig_18 = sys.argv[4]
    contig_36 = sys.argv[5]
    output_fasta = sys.argv[6]

    circularize_assembly(assembly_fasta, bridge_fasta, contig_8, contig_18, contig_36, output_fasta)