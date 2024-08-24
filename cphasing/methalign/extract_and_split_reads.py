#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2024-05-27 15:30

import argparse
import pysam
import re
import os
import collections

pysam.set_verbosity(0)

comp_table = str.maketrans('ATCG', 'TAGC')


def rev_com(seq):

    seq_upper = seq.upper()
    assert 'N' not in seq_upper and '-' not in seq_upper
    
    return seq_upper.translate(comp_table)[::-1]


def get_mm(aln, query_start, c_index_list, query_aln_len, split_len, read_seq):
    
    mm_tag = aln.get_tag('MM')
    ml_tag = aln.get_tag('ML')

    mm_tag_prefix = mm_tag.split(',', 1)[0]
    nsubreads = query_aln_len // split_len
    sub_mm_list = [[] for _ in range(nsubreads)]
    sub_ml_list = [[] for _ in range(nsubreads)]
    
    index_sum = 0
    for n, i in enumerate((int(shift) for shift in mm_tag.rstrip(';').split(',')[1:])):
        index_sum += i
        subread_index = (c_index_list[n+index_sum] - query_start) // split_len
        if subread_index < 0 or subread_index > nsubreads - 1:
            continue
        if sub_mm_list[subread_index]:
            sub_mm_list[subread_index].append(str(i))
        else:
            sub_mm_list[subread_index].append(str(read_seq[subread_index*split_len+query_start: c_index_list[n+index_sum]].count('C')))
        sub_ml_list[subread_index].append(str(ml_tag[n]))
    
    sub_mm_tags = ['{},{};'.format(mm_tag_prefix, ','.join(sub_mm)) for sub_mm in sub_mm_list]
    sub_ml_tags = [','.join(sub_ml) for sub_ml in sub_ml_list]
    
    return sub_mm_tags, sub_ml_tags


def construct_subreads(aln, split_len, header):

    if aln.is_forward:
        read_seq = aln.query_sequence
        seq = aln.query_alignment_sequence
        qual = ''.join([chr(q+33) for q in aln.query_alignment_qualities])
        if aln.cigartuples[0][0] == 5:
            hard_clip = aln.cigartuples[0][1]
        else:
            hard_clip = 0
        query_start = aln.query_alignment_start
        aligned_pairs = [(q, r) for q, r in aln.get_aligned_pairs() if q is not None and r is not None]
    else:
        read_seq = rev_com(aln.query_sequence)
        seq = rev_com(aln.query_alignment_sequence)
        qual = ''.join([chr(q+33) for q in aln.query_alignment_qualities[::-1]])
        if aln.cigartuples[-1][0] == 5:
            hard_clip = aln.cigartuples[-1][1]
        else:
            hard_clip = 0
        read_len = aln.infer_read_length()
        query_start = read_len - aln.query_alignment_end + hard_clip
        aligned_pairs = [(read_len - q - 1, r) for q, r in aln.get_aligned_pairs()[::-1] if q is not None and r is not None]
    
    ref_range_list = [[-1, -1] for _ in range(aln.query_alignment_length // split_len)]
    
    for q, r in aligned_pairs:
        subread_index = (q - query_start) // split_len
        if subread_index >= 0 and subread_index < len(ref_range_list):
            if ref_range_list[subread_index][0] == -1 or r < ref_range_list[subread_index][0]:
                ref_range_list[subread_index][0] = r
            if ref_range_list[subread_index][1] == -1 or r > ref_range_list[subread_index][1]:
                ref_range_list[subread_index][1] = r

    matches = re.finditer(r'C', read_seq)
    c_index_list = [match.start() for match in matches]
    sub_mm_tags, sub_ml_tags = get_mm(aln, query_start, c_index_list, aln.query_alignment_length, split_len, read_seq)
    subreads = []
    
    for n, start in enumerate(range(0, aln.query_alignment_length, split_len)):
        if start + split_len <= aln.query_alignment_length and ref_range_list[n][0] != -1 and ref_range_list[n][1] != -1:
            sub_seq = seq[start:start+split_len]
            sub_qual = qual[start:start+split_len]
            orient = '+' if aln.is_forward else '-'
            subread_name = '{}_{}_{}_{}_{}_{}_{}'.format(
                    aln.query_name, split_len, n, aln.reference_name, ref_range_list[n][0]+1, ref_range_list[n][1]+1, orient)
            if sub_ml_tags[n]:
                subread = pysam.AlignedSegment.fromstring('{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tRG:Z:{}\tMM:Z:{}\tML:B:C,{}'.format(
                    subread_name, sub_seq, sub_qual, aln.get_tag('RG'), sub_mm_tags[n], sub_ml_tags[n]), header)
            else:
                subread = pysam.AlignedSegment.fromstring('{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tRG:Z:{}'.format(
                    subread_name, sub_seq, sub_qual, aln.get_tag('RG')), header)
            subreads.append(subread)
    
    return subreads


def get_header(bam_reads):
    
    for n, bam in enumerate(bam_reads):
        with pysam.AlignmentFile(bam, mode='rb', check_sq=False) as f:
            if n == 0:
                header_dict = f.header.as_dict()
            else:
                header_dict['RG'].extend(f.header.as_dict()['RG'])
    
    return pysam.AlignmentHeader.from_dict(header_dict)


def parse_bam(bam_alignments, bam_reads, split_len, min_len, threads):
    
    out_bam = 'split_{}_{}_'.format(min_len, split_len) + os.path.basename(bam_alignments)
    with pysam.AlignmentFile(bam_alignments, mode='rb', format_options=[b'filter=!flag.unmap'], check_sq=False, threads=threads) as fin:
        header = get_header(bam_reads)
        with pysam.AlignmentFile(out_bam, mode='wb', header=header, threads=threads) as fout:
            for aln in fin:
                # only primary alignments with a length >= min_len is used for subread splitting
                if aln.query_alignment_length < min_len or aln.is_supplementary or aln.is_secondary:
                    continue
                subreads = construct_subreads(aln, split_len, header)
                for subread in subreads:
                    fout.write(subread)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bam_alignments', help='input UL read alignments in BAM format')
    parser.add_argument('bam_reads', nargs='+', help='bam file containing raw UL reads (for header reconstruction)')
    parser.add_argument('--min_len', type=int, default=100000, help='minimum length of reads for splitting subreads, default: %(default)s bp')
    parser.add_argument('--split_len', type=int, default=250, help='length of split subreads, default: %(default)s bp')
    parser.add_argument('--threads', type=int, default=8, help='threads for parsing BAM file, default: %(default)s')
    args = parser.parse_args()

    parse_bam(args.bam_alignments, args.bam_reads, args.split_len, args.min_len, args.threads)


if __name__ == '__main__':
    main()

