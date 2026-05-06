#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse
import logging
import os
import os.path as op
import sys


import sys
from collections import defaultdict

def extract_prefix(contig_name):

    return contig_name.split('.')[0]

def main(input_agp, output_agp):
    group_prefix_lengths = defaultdict(lambda: defaultdict(int))
    skip_groups = set()
    with open(input_agp, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            
            parts = line.split('\t')
            if len(parts) >= 8 and parts[4] == 'W':
                group_id = parts[0]
                contig_name = parts[5]
                if group_id == contig_name:
                    skip_groups.add(group_id)
                    continue
                length = int(parts[7])
                
                prefix = extract_prefix(contig_name)
                group_prefix_lengths[group_id][prefix] += length

    group_to_new_name = {}
    prefix_usage_count = defaultdict(int)
    
    for group_id, prefix_counts in group_prefix_lengths.items():
        best_prefix = sorted(prefix_counts.items(), key=lambda x: x[1], reverse=True)[0][0]
        
        prefix_usage_count[best_prefix] += 1
        if prefix_usage_count[best_prefix] == 1:
            new_name = best_prefix
        else:
            new_name = f"{best_prefix}_{prefix_usage_count[best_prefix]}"
            
        group_to_new_name[group_id] = new_name
        
        print(f"Group: {group_id} -> Assigned: {new_name} (Majority: {best_prefix})")

    with open(input_agp, 'r') as fin, open(output_agp, 'w') as fout:
        for line in fin:
            if not line.strip() or line.startswith("#"):
                fout.write(line)
                continue
            
            parts = line.rstrip('\n').split('\t')
            old_group = parts[0]
            
            if old_group in group_to_new_name:
                parts[0] = group_to_new_name[old_group]
                
            fout.write('\t'.join(parts) + '\n')

    print(f"\nDone! Renamed AGP has been saved to: {output_agp}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rename_agp_by_majority.py <input.agp> <output.agp>")
        sys.exit(1)
        
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    main(input_file, output_file)