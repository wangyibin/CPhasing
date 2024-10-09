#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2024-02-05 17:04


import argparse
import sys
import os
import pysam
import re
from array import array
import collections
import numpy as np
from portion import closed


# from line_profiler import profile


pysam.set_verbosity(0)


comp_table = str.maketrans('ATCG', 'TAGC')


def parse_fasta(fasta):
    
    fa_dict = dict()
    with open(fasta) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                ID = line.split()[0][1:]
                fa_dict[ID] = []
            else:
                fa_dict[ID].append(line.strip().upper())

    for ID, seq_list in fa_dict.items():
        fa_dict[ID] = ''.join(seq_list)

    return fa_dict


def parse_bed(bed, cov_cutoff=75):
    
    hifi_methylation_dict = collections.defaultdict(set)
    with open(bed) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.split()
            ctg, pos, cov = cols[0], int(cols[1]), float(cols[8])
            if cov >= cov_cutoff:
                hifi_methylation_dict[ctg].add(pos)

    return hifi_methylation_dict

# @profile
def get_5mc_sites_from_read(read_name, read_seq, MM_tag, ML_tag, prob_cutoff, pattern=r'C'):
   
    matches = re.finditer(pattern, read_seq)
    ont_pos_index = {match.start(): index for index, match in enumerate(matches)}
    ont_methylation_array = array('h', [0] * len(ont_pos_index))
    index_sum = 0
    for n, i in enumerate((int(shift) for shift in MM_tag.rstrip(';').split(',')[1:])):
        index_sum += i
        if ML_tag[n] >= prob_cutoff:
            ont_methylation_array[n+index_sum] = 1
            
    return ont_pos_index, ont_methylation_array

# @profile
def get_query_alignment_termini(aln):

    if aln.is_forward:
        # '5' means hard-clipping
        if aln.cigartuples[0][0] == 5:
            hard_clip = aln.cigartuples[0][1]
        else:
            hard_clip = 0
        return aln.query_alignment_start + hard_clip, aln.query_alignment_end + hard_clip
    else:
        if aln.cigartuples[-1][0] == 5:
            hard_clip = aln.cigartuples[-1][1]
        else:
            hard_clip = 0
        read_len = aln.infer_read_length()
        return read_len - aln.query_alignment_end + hard_clip, read_len - aln.query_alignment_start + hard_clip
    
# @profile    
def condense_cigar(cigartuples):
    
    lsoft, rsoft, match = 0, 0, 0
    insertion, deletion = 0, 0
    for operation, length in cigartuples:
        if operation == 0:
            match += length
        elif operation == 1:
            insertion += length
        elif operation == 2:
            deletion += length
        elif operation not in {4, 5}:
            raise Exception('Unknown cigar operation: {}'.format(operation))
    if cigartuples[0][0] in {4, 5}:
        lsoft += cigartuples[0][1]
    if cigartuples[-1][0] in {4, 5}:
        rsoft += cigartuples[-1][1]
    lsoft_str = '{}S'.format(lsoft) if lsoft else ''
    rsoft_str = '{}S'.format(rsoft) if rsoft else ''
    if insertion > deletion:
        match_str = '{}M'.format(match + deletion)
    else:
        match_str = '{}M'.format(match + insertion)
    indel = deletion - insertion
    if indel:
        if indel > 0:
            indel_str = '{}D'.format(indel)
        else:
            indel_str = '{}I'.format(-indel)
    else:
        indel_str = ''
    return '{}{}{}{}'.format(lsoft_str, match_str, indel_str, rsoft_str)


# @profile
def reconstruct_SA_tag(aln):
    
    strand = '+' if aln.is_forward else '-'
    cigarstring = condense_cigar(aln.cigartuples)
    SA_tag = '{},{},{},{},{},{};'.format(
            aln.reference_name, aln.reference_start + 1, 
            strand, cigarstring,
            aln.mapping_quality, 
            aln.get_tag('NM'))

    return SA_tag


def parse_bam(bam, fa_dict, hifi_methylation_dict, penalty, prob_cutoff, designate_mapq, recalculate_all, threads):
    
    print('Read\tRef\tquery_aln_len\tflag\tMAPQ\tPos\tAlignment_score\tInconsistent_5mC\tAdjusted_score')
    
    primary_flags = {0, 16}
    cached_aln = []

    # @profile
    def recalculate_score(aln, flg, ont_pos_index, ont_methylation_array):
        inconsistent_5mC = 0
        ref_name, is_forward = aln.reference_name, aln.is_forward
        hifi_methylation_set = hifi_methylation_dict[ref_name]
        fa_seq = fa_dict[ref_name]
        for read_aln_pos, ref_aln_pos in aln.get_aligned_pairs():
            raw_read_aln_pos = read_aln_pos
            # for reverse strand, recalculate read_aln_pos based on the contig length
            if not is_forward and read_aln_pos is not None:
                read_aln_pos = aln.infer_query_length() - (read_aln_pos + 1)
            # not C in ONT read
            if read_aln_pos not in ont_pos_index:
                continue
            
            is_methylated = ont_methylation_array[ont_pos_index[read_aln_pos]]
            
            # forward strand
            if is_forward:
                # mismatch
                if ref_aln_pos is None or fa_seq[ref_aln_pos] != 'C':
                    continue
                if not ((ref_aln_pos in hifi_methylation_set and is_methylated) or 
                        (ref_aln_pos not in hifi_methylation_set and not is_methylated)):
                    inconsistent_5mC += 1
            # reverse strand
            else:
                # mismatch
                if ref_aln_pos is None or fa_seq[ref_aln_pos] != 'G':
                    continue
                if not ((ref_aln_pos - 1 in hifi_methylation_set and is_methylated) or
                        (ref_aln_pos - 1 not in hifi_methylation_set and not is_methylated)):
                    inconsistent_5mC += 1
        
        # calculate adjusted score
        score = aln.get_tag('AS')
        adjusted_score = score - inconsistent_5mC * penalty
        aln.set_tag('AS', adjusted_score)
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
            last_read, ref_name, aln.query_alignment_length, flg,
            aln.mapping_quality, aln.pos, score, inconsistent_5mC, adjusted_score))
    
    # @profile
    def analyze_cached_aln(fout, last_read, alignment, flag):
        
        if len(cached_aln) > 1:
            
            supplementary_alns = []
            non_supplementary_alns = []
            aln_list = []
            any_primary, any_supplementary, any_secondary = False, False, False
            
            for aln_flg in cached_aln:
                aln, flg = aln_flg
                if flg in primary_flags or aln.is_supplementary:
                    aln_list.append([aln_flg])
                    if flg in primary_flags:
                        any_primary = True
                        # get necessary information from the primary alignment
                        read_length = aln.infer_query_length()
                        read_seq = aln.query_sequence if aln.is_forward else aln.query_sequence.translate(comp_table)[::-1]
                        if aln.has_tag('MM'):
                            assert aln.has_tag('ML')
                            mm_tag, ml_tag = aln.get_tag('MM'), aln.get_tag('ML')
                        else:
                            mm_tag, ml_tag = '', ''
                        ont_pos_index, ont_methylation_array = get_5mc_sites_from_read(last_read, read_seq, mm_tag, ml_tag, prob_cutoff)
                    if aln.is_supplementary:
                        any_supplementary = True
                else:
                    assert aln.is_secondary
                    any_secondary = True
            
            # skip the read if there is no primary alignment
            if any_primary:
                if any_supplementary and any_secondary:
                    ps_qry_range_list = []
                    for alns in aln_list:
                        ps_qry_start, ps_qry_end = get_query_alignment_termini(alns[0][0])
                        ps_qry_range_list.append((ps_qry_start, ps_qry_end))
                
                if any_secondary:
                    for aln, flg in cached_aln:
                        # find best assignment of secondary alignments
                        if any_supplementary and aln.is_secondary:
                            ovl_len_list = []
                            # recalculate score for primary and supplementary alignments (always necessary)
                            recalculate_score(aln, flg, ont_pos_index, ont_methylation_array)
                            for n, (ps_qry_start, ps_qry_end) in enumerate(ps_qry_range_list):
                                query_start, query_end = get_query_alignment_termini(aln)
                                overlap = closed(query_start + 1, query_end) & closed(ps_qry_start + 1, ps_qry_end)
                                ovl_len = overlap.upper - overlap.lower + 1 if overlap else 0
                                ovl_len_list.append((n, ovl_len))

                            ovl_len_list.sort(key=lambda x: x[1], reverse=True)
                            best_ps_index, best_ps_ovl = ovl_len_list[0]
                            aln_list[best_ps_index].append((aln, flg))
                        elif aln.is_secondary:
                            recalculate_score(aln, flg, ont_pos_index, ont_methylation_array)
                            aln_list[0].append((aln, flg))
                    # recalculate score for primary and supplementary alignments
                    # if there are no corresponding secondary alignments, recalculate the score of the primary/supplementary
                    # alignment is not necessary, unless --recalculate_all is set
                    for alns in aln_list:
                        if len(alns) > 1 or recalculate_all:
                            ps_aln, ps_flg = alns[0]
                            recalculate_score(ps_aln, ps_flg, ont_pos_index, ont_methylation_array)

                    # if any_secondary:
                    SA_tag_list = []
                    output_aln_list = []
                    any_change = False
                    aln_list.sort(key=lambda x: x[0][1])
                    for alns in aln_list:
                        original_best_aln = alns[0][0]
                        # no corresponding secondary alignments for this primary/supplementary alignment
                        if len(alns) == 1:
                            output_aln_list.append([original_best_aln])
                            if any_supplementary:
                                SA_tag_list.append(reconstruct_SA_tag(original_best_aln))
                            continue
                        # find the best alignment based on adjusted scores
                        # if a primary/supplementary alignment has the same adjusted score as a secondary counterpart, 
                        # choose the primary/supplementary alignment as the best one
                        alns.sort(key=lambda x: x[0].get_tag('AS'), reverse=True)
                        best_aln = alns[0][0]
                        # the original best alignment is still the best after recalculating the scores
                        if best_aln == original_best_aln:
                            output_aln_list.append([aln for aln, _ in alns])
                            if any_supplementary:
                                SA_tag_list.append(reconstruct_SA_tag(best_aln))
                            continue
                        any_change = True
                        # check
                        for n, (aln, flg) in enumerate(alns):
                            if (flg in primary_flags or aln.is_supplementary):
                                ps_aln = aln
                                if aln.get_tag('AS') == best_aln.get_tag('AS'):
                                    assert n == 0
                        # TODO: reconstruct query_squence and query_qualities according to raw reads 
                        assert best_aln.is_secondary
                        if best_aln.is_forward:
                            if ps_aln.flag in primary_flags:
                                best_aln.flag = 0
                            else:
                                assert ps_aln.is_supplementary
                                best_aln.flag = 2048
                        else:
                            if ps_aln.flag in primary_flags:
                                best_aln.flag = 16
                            else:
                                assert ps_aln.is_supplementary
                                best_aln.flag = 2064
                        if designate_mapq != -1:
                            best_aln.mapping_quality = designate_mapq
                        
                        output_aln_list.append([best_aln])
                        
                        if any_supplementary:
                            SA_tag_list.append(reconstruct_SA_tag(best_aln))
                        
                        # other alignments -> secondary
                        for aln, flg in alns[1:]:
                            if flg in primary_flags or aln.is_supplementary:
                                # modify some fields for the original primary/supplementary alignments
                                aln.query_sequence = None
                                aln.query_qualities = None
                                if aln.is_forward:
                                    aln.flag = 256
                                else:
                                    aln.flag = 272
                                aln.mapping_quality = 0
                                aln.set_tag('SA', None)
                            output_aln_list[-1].append(aln)
                    
                    # reconstruct the SA tags for new primary/supplementary alignments
                    if any_supplementary and any_change:
                        for n, alns in enumerate(output_aln_list):
                            ps_aln = alns[0]
                            SA_tag = ''.join([tag for i, tag in enumerate(SA_tag_list) if i != n])
                            ps_aln.set_tag('SA', SA_tag, value_type='Z')
                            for aln in alns:
                                fout.write(aln)
                    else:
                        for alns in output_aln_list:
                            for aln in alns:
                                fout.write(aln)
                else:
                    cached_aln.sort(key=lambda x: x[1])
                    for aln, flg in cached_aln:
                        if recalculate_all:
                            recalculate_score(aln, flg, ont_pos_index, ont_methylation_array)
                        fout.write(aln)
        else:
            aln, flg = cached_aln[0]
            if flg in primary_flags:
                if recalculate_all:
                    read_seq = aln.query_sequence if aln.is_forward else ''.join(aln.query_sequence.translate(comp_table)[::-1])
                    if aln.has_tag('MM'):
                        assert aln.has_tag('ML')
                        mm_tag, ml_tag = aln.get_tag('MM'), aln.get_tag('ML')
                    else:
                        mm_tag, ml_tag = '', ''
                    ont_pos_index, ont_methylation_array = get_5mc_sites_from_read(last_read, read_seq, mm_tag, ml_tag, prob_cutoff)
                    recalculate_score(aln, flg, ont_pos_index, ont_methylation_array)
                else:
                    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\twhatever'.format(
                        last_read, aln.reference_name, aln.query_alignment_length, flg, aln.mapping_quality, aln.pos, aln.get_tag('AS')))
                fout.write(aln)
        last_read = read_name
        cached_aln.clear()
        cached_aln.append([alignment, flag])
    
        return last_read

    base_bam = os.path.basename(bam)
    out_bam = '{}.penalty{}.{}'.format(os.path.splitext(base_bam)[0], penalty, 'methly_filtered.bam')
    
    format_options = [b'filter=!flag.unmap']
    
    with pysam.AlignmentFile(bam, format_options=format_options, threads=threads) as fin:
        with pysam.AlignmentFile(out_bam, 'wb', template=fin, threads=threads) as fout:
            last_read = ''
            for alignment in fin:
                read_name, flag = alignment.query_name, alignment.flag
                if read_name == last_read:
                    cached_aln.append([alignment, flag])
                elif last_read:
                    last_read = analyze_cached_aln(fout, last_read, alignment, flag)
                else:
                    last_read = read_name
                    cached_aln.append([alignment, flag])
            analyze_cached_aln(fout, last_read, alignment, flag)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta file')
    parser.add_argument('bed', help='bed file output by pb-CpG-tools')
    parser.add_argument('bam', help='input bam file of ONT read mapping results')
    parser.add_argument('--penalty', type=int, default=2, help='penalty for inconsistent 5mC, default: %(default)s')
    parser.add_argument('--prob_cutoff', type=int, default=128, help='probability cutoff, default: %(default)s')
    parser.add_argument('--designate_mapq', type=int, default=60, help='designate MAPQ for the best alignments. Set the value to -1 to disable the function, default: %(default)s')
    parser.add_argument('--recalculate_all', default=False, action='store_true', help='recalculate scores for all alignments although there may not be any secondary alignments, default: %(default)s')
    parser.add_argument('--threads', type=int, default=8, help='number of threads for bam reading, default: %(default)s')
    args = parser.parse_args()
    
    if not -1 <= args.designate_mapq <= 60:
        raise Exception('--designate_mapq should be within the range [-1, 60]')
    
    fa_dict = parse_fasta(args.fasta)
    hifi_methylation_dict = parse_bed(args.bed)
    parse_bam(args.bam, fa_dict, hifi_methylation_dict, args.penalty, args.prob_cutoff, args.designate_mapq, args.recalculate_all, args.threads)


if __name__ == '__main__':
    main()

