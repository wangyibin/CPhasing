#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2022-04-03 20:27

import os
import sys
import argparse
import collections
import glob

from collections import defaultdict
from pathlib import Path
from tempfile import TemporaryDirectory
from cphasing.agp import import_agp

def parse_fasta(fasta):
    
    fa_len_dict = dict()
    ignore = False
    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                if 'chimeric' in line:
                    ignore = True
                    continue
                else:
                    ignore = False
                    fa_len_dict[ID] = 0

            elif not ignore:
                fa_len_dict[ID] += len(line.strip())
            
    return fa_len_dict


def parse_groups(groups, fa_len_dict):

    n_groups = len(groups)
    anchored_len_dict = collections.defaultdict(int)
    inter_homo_err = 0
    inter_nonhomo_err = 0
    
    chr4_8_excluded_anchored_len = 0
    chr4_8_excluded_inter_homo_err = 0
    chr4_8_excluded_inter_nonhomo_err = 0

    largest_group_dict = collections.defaultdict(int)
    inter_homo_contigs = []
    for group in groups.values.tolist():
        
        source_chr_len_dict = collections.defaultdict(int)
        chr4_8_excluded_group_len = 0

        nctgs = len(group)

        most_group_length_db = defaultdict(int)
        ctg_info_list = []
        for ctg in group:
            if 'collapsed' in ctg or 'chimeric' in ctg:
                continue 
            length = fa_len_dict[ctg]
            source_chr = ctg.split('.')[0]
            
            most_group_length_db[source_chr] += length 
            ctg_info_list.append((ctg, source_chr, length))

        most_group = max(most_group_length_db, key=lambda x: most_group_length_db[x])


        # one group should have at least two contigs inside
        if nctgs < 2:
            print('group file {} is skipped because of {} contig inside'.format(
                group, nctgs), file=sys.stderr)
            continue

        source_chr_len_dict.clear()

        for ctg in group:
            # collapsed and chimeric contigs are skipped
          
            if 'chimeric' in ctg:
                    continue
          
            length = fa_len_dict[ctg]
            if 'collapsed' in ctg:
                ctg1 = ctg.split(":")[0]
                ctg2 = ctg.split("|")[1].split(":")[0]
                source_chr1 = ctg1.split('.')[0]
                source_chr2 = ctg2.split('.')[0]

                if source_chr1 == most_group:
                    source_chr = source_chr1
                elif source_chr2 == most_group:
                    source_chr = source_chr2 
               
                anchored_len_dict[source_chr] += length

            else:
                source_chr = ctg.split('.')[0]
                anchored_len_dict[source_chr] += length
            
            source_chr_len_dict[source_chr] += length
           
            if source_chr.split('_')[0] not in {'Chr4', 'Chr8'}:
                chr4_8_excluded_group_len += length
    
        source_chr_len_list = list(source_chr_len_dict.items())
        if source_chr_len_list:
            source_chr_len_list.sort(key=lambda x: x[1])
            dominant_chr, dominant_chr_len = source_chr_len_list[-1]
        else:
            continue

        for source_chr, length in source_chr_len_list:
            if length > largest_group_dict[source_chr]:
                largest_group_dict[source_chr] = length

        for ctg, source_chr, length in ctg_info_list:
            if source_chr == dominant_chr:
                continue
            if source_chr[0] == dominant_chr[0]:
                inter_homo_err += length
                inter_homo_contigs.append(
                    (ctg, source_chr, dominant_chr, length)
                )
            else:
                inter_nonhomo_err += length
    

    return n_groups, anchored_len_dict, inter_homo_err, inter_nonhomo_err, largest_group_dict, inter_homo_contigs


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('agp', help='input agp file')
    # parser.add_argument('--real', help="real chrom and contig pair", default=None)
    args = parser.parse_args()

    agp = str(Path(args.agp).absolute())


    agp_df, _ = import_agp(agp)
    agp_df.reset_index(inplace=True)
    contig_sizes = agp_df[['id', 'tig_end']]
    contig_sizes = contig_sizes.set_index('id')

    agp_df = agp_df[agp_df['chrom'] != agp_df['id']]
    # 
    agp_df = agp_df[~agp_df['chrom'].str.contains('ctg')]
    agp_df['chrom'] = agp_df['chrom'].astype(str)
    cluster_df = agp_df.groupby('chrom')['id'].apply(lambda x: list(x))


    contigsizes = contig_sizes.to_dict()['tig_end']
    fa_len_dict = {}
    for contig, length in contigsizes.items():
        if 'chimeric' in contig:
            continue 

        fa_len_dict[contig] = length
    
    fa_len_dict = dict(map(lambda x: (x[0], int(x[1])), fa_len_dict.items()))

    

    total_ctg_len = 0
    for ctg in fa_len_dict:
        if 'collapsed' in ctg:
            total_ctg_len += fa_len_dict[ctg] * 2
        else:
            total_ctg_len += fa_len_dict[ctg]
    
    n_groups, anchored_len_dict, inter_homo_err, inter_nonhomo_err, largest_group_dict, inter_homo_contigs = parse_groups(cluster_df, fa_len_dict)

    anchored_len = sum(anchored_len_dict.values())

    try:
        contiguity = sum([d_chr_len/anchored_len_dict[d_chr] for d_chr, d_chr_len in largest_group_dict.items()]) / len(anchored_len_dict)
    except:
        contiguity = 0
    
    
    inter_homo_error_rate = inter_homo_err/anchored_len if anchored_len else 0 

    inter_nonhomo_error_rate = inter_nonhomo_err/anchored_len  if anchored_len else 0
    print('Contiguity\t{}'.format(contiguity))
    print('Inter_homo_error_rate(%)\t{}'.format(inter_homo_error_rate))
    print('Inter_nonhomo_error_rate(%)\t{}'.format(inter_nonhomo_error_rate))
    print('Ngroups\t{}'.format(n_groups))
    print('Anchoring rate(%)\t{}'.format(anchored_len/total_ctg_len*100))

    out_tsv = os.path.splitext(agp)[0] + ".inter_homo_contigs.tsv"
    with open(out_tsv, "w") as fh:
        print("contig\tcurrent_chr\tcorrect_chr\tlength", file=fh)
        for ctg, cur_chr, cor_chr, length in inter_homo_contigs:
            print(f"{ctg}\t{cur_chr}\t{cor_chr}\t{length}", file=fh)
    print(f"Inter-homo error contig list written to {out_tsv}", file=sys.stderr)
if __name__ == '__main__':
    main()
