#!/usr/bin/env python3
import gzip
import sys

input_file = "wRi_Merrill_23_annotated.embl.gz"
output_file = "wRi_Merrill_23_annotated_fixed.embl.gz"

fixed_count = 0

with gzip.open(input_file, 'rt') as f_in, gzip.open(output_file, 'wt') as f_out:
    for line_num, line in enumerate(f_in, 1):
        # Fix qualifier lines that are missing FT prefix
        # These start with exactly 19 spaces followed by /
        if line.startswith('                   /'):
            line = 'FT' + line[2:]
            fixed_count += 1
            if fixed_count <= 10:  # Print first 10 fixes
                print(f"Fixed line {line_num}", file=sys.stderr)
        
        f_out.write(line)

print(f"\nTotal lines fixed: {fixed_count}", file=sys.stderr)
print("Done!", file=sys.stderr)
