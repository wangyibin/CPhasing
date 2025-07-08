#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2024-02-05 17:04


import argparse
import logging
import sys
import os
import pysam
import re
from array import array
import collections
import numpy as np
import pandas as pd 
import polars as pl
from portion import closed

from collections import deque, defaultdict

logger = logging.getLogger(__name__)

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

def parse_bed_modkit(bed, cov_cutoff=75):
    
    hifi_methylation_dict = defaultdict(set)
    with open(bed) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.split()
            ctg, pos, cov = cols[0], int(cols[1]), float(cols[10])
            if cov >= cov_cutoff:
                hifi_methylation_dict[ctg].add(pos)

    return hifi_methylation_dict


def parse_bedgraph(bedgraph, cov_cutoff=75, threads=4):
    os.environ["POLARS_MAX_THREADS"] = str(threads)

    df = pl.read_csv(
        bedgraph,
        separator='\t',
        has_header=False,
        new_columns=['ctg', 'start', 'end', 'cov'],
        comment_prefix='#',
        dtypes={'ctg': pl.Utf8, 'start': pl.Int64, 'end': pl.Int64, 'cov': pl.Float64}
    ).select(
        pl.col('ctg'),
        pl.col('start').cast(pl.Int64),
        pl.col('cov').cast(pl.Float64)
    ).filter(pl.col('cov') >= cov_cutoff)
   
    gb = df.group_by('ctg').agg(
        pl.col('start').unique().alias('starts')
    )
    hifi_methylation_dict = defaultdict(set,
        {ctg: set(starts) for ctg, starts in zip(gb['ctg'], gb['starts'])}
    )

    return hifi_methylation_dict


def get_5mc_sites_from_read(read_seq, MM_tag, ML_tag, prob_cutoff, pattern=re.compile(r'C')):
    if read_seq is None or not read_seq:
        return {}, array('h', [])
    matches = pattern.finditer(read_seq)
    ont_pos_index = {match.start(): index for index, match in enumerate(matches)}
    ont_methylation_array = array('h', [0] * len(ont_pos_index))
    index_sum = 0

    # _MM_tags = MM_tag.rstrip(';').split(",") 
    # idx = [idx if 'C+m' in val else None for idx, val in enumerate(_MM_tags)]
    # idx = list(filter(lambda x: x is not None, idx ))
    # if not idx:
    #     return ont_pos_index, ont_methylation_array
    # else:
    #     idx = idx[0]

    # MM_tags = MM_tag.rstrip(';').split(';')
    # MM_tags = list(filter(lambda x: 'C+m' in x.split(";")[0], MM_tags))
    
    # if len(MM_tags) == 0:
    #     return ont_pos_index, ont_methylation_array
    # MM_tags = MM_tags[0].split(",")[1:]
    # MM_tags_count = len(MM_tags)
    # MM_tags = list(map(int, MM_tags))
    # ML_tag = ML_tag[idx: idx + MM_tags_count]

    parts = MM_tag.rstrip(';').split(';')
    offset = 0
    shifts = None
    for part in parts:
        items = part.split(',')
        if 'C+m' in items[0]:
            shifts = list(map(int, items[1:]))
            break
        offset += len(items)
    if not shifts:
        return ont_pos_index, ont_methylation_array
    MM_tags = shifts
   
    ML_tag = ML_tag[offset: offset + len(shifts)]

    for n, i in enumerate((int(shift) for shift in MM_tags)):
        index_sum += i
        if ML_tag[n] >= prob_cutoff:
            try:
                ont_methylation_array[n+index_sum] = 1
            except IndexError:
                continue
            
    return ont_pos_index, ont_methylation_array

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
    
    # print('Read\tRef\tquery_aln_len\tflag\tMAPQ\tPos\tAlignment_score\tInconsistent_5mC\tAdjusted_score')
    
    primary_flags = {0, 16}
    def recalculate_score(aln, flg, ont_pos_index, ont_methylation_array):
        inconsistent_5mC = 0
        ref_name, is_forward = aln.reference_name, aln.is_forward
        hifi_set = hifi_methylation_dict[ref_name]
        fa_seq = fa_dict[ref_name]
        aligned_pairs = aln.get_aligned_pairs()
        # build a map from read-side positions to reference positions,
        # applying reverse-strand adjustment once
        # read_to_ref = {}
        # read_len = aln.infer_query_length()

        # read_to_ref = [-1] * read_len
        # for rpos, refpos in aligned_pairs:
        #     if rpos is None or refpos is None:
        #         continue
        #     if not is_forward:
        #         rpos = read_len - (rpos + 1)
        #     read_to_ref[rpos] = refpos

        # # now iterate only over C positions in the read (ont_pos_index keys)
        # for rpos, idx in ont_pos_index.items():
        #     if rpos < 0 or rpos >= read_len:
        #         continue
        #     refpos = read_to_ref[rpos]
        #     if refpos < 0:
        #         continue
        #     is_meth = bool(ont_methylation_array[idx])

        #     if is_forward:
        #         if fa_seq[refpos] != 'C':
        #             continue
        #         hit = refpos in hifi_set
        #     else:
        #         if fa_seq[refpos] != 'G':
        #             continue
        #         hit = (refpos - 1) in hifi_set

        #     if hit != is_meth:
        #         inconsistent_5mC += 1
        read_len = aln.infer_query_length()
        _idx = ont_pos_index  # local ref
        _meth = ont_methylation_array
        _hifi = hifi_set
        seq = fa_seq
        inc = 0
        for rpos, refpos in aligned_pairs:
            if rpos is None or refpos is None:
                continue
            # adjust for reverse strand
            if not is_forward:
                rpos = read_len - (rpos + 1)
            # only care positions with C in read
            idx = _idx.get(rpos)            
            if idx is None:
                    continue
            is_meth = _meth[idx] != 0
            if is_forward:
                if seq[refpos] != 'C':
                    continue
                hit = (refpos in _hifi)
            else:
                if seq[refpos] != 'G':
                    continue
                hit = ((refpos - 1) in _hifi)
            if hit != is_meth:
                inc += 1
        inconsistent_5mC = inc
        # update tags
        orig_score = aln.get_tag('AS')
        new_score = orig_score - inconsistent_5mC * penalty
        aln.set_tag('MA', inconsistent_5mC)
        aln.set_tag('AS', new_score)
  
    def analyze_cached_aln(cached_aln, flag):
        output_alns = []
        if len(cached_aln) > 1:
            aln_list = []
            any_primary, any_supplementary, any_secondary = False, False, False
            
            for aln_flg in cached_aln:
                aln, flg = aln_flg
                if flg in primary_flags or aln.is_supplementary:
                    aln_list.append([aln_flg])
                    if flg in primary_flags:
                        any_primary = True
                        # get necessary information from the primary alignment
                        # read_length = aln.infer_query_length()
                        read_seq = aln.query_sequence if aln.is_forward else aln.query_sequence.translate(comp_table)[::-1]
                        if aln.has_tag('MM'):
                            assert aln.has_tag('ML')
                            mm_tag, ml_tag = aln.get_tag('MM'), aln.get_tag('ML')
                        else:
                            mm_tag, ml_tag = '', ''
                        ont_pos_index, ont_methylation_array = get_5mc_sites_from_read(read_seq, mm_tag, ml_tag, prob_cutoff)
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
                    sec_infos = []
                    for aln, flg in cached_aln:
                        if aln.is_secondary:
                            # recalculate_score(aln, flg, ont_pos_index, ont_methylation_array)
                            qs, qe = get_query_alignment_termini(aln)
                            sec_infos.append((aln, flg, qs, qe))

                    for aln, flg, qs, qe in sec_infos:
                        if any_supplementary:
                            ovl_len_list = []
                            for i, (ps_start, ps_end) in enumerate(ps_qry_range_list):
                                overlap = closed(qs+1, qe) & closed(ps_start+1, ps_end)
                                ovl = overlap.upper - overlap.lower + 1 if overlap else 0
                                ovl_len_list.append((i, ovl))
                            best_group = max(ovl_len_list, key=lambda x: x[1])[0]
                            aln_list[best_group].append((aln, flg))
                        else:
                            aln_list[0].append((aln, flg))

                    # recalculate score for primary and supplementary alignments
                    # if there are no corresponding secondary alignments, recalculate the score of the primary/supplementary
                    # alignment is not necessary, unless --recalculate_all is set
                    for alns in aln_list:
                        if len(alns) > 1 or recalculate_all:
                            ps_aln, ps_flg = alns[0]
                            for aln, flg in alns:
                                recalculate_score(aln, flg, ont_pos_index, ont_methylation_array)

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
                         
                            if designate_mapq != -1 and best_aln.mapping_quality < designate_mapq:
                                if len(alns) >= 2:
                                    second_aln = alns[1][0]
                                    if second_aln.get_tag('AS') == best_aln.get_tag('AS'):
                                        continue
                                best_aln.mapping_quality = designate_mapq
                                best_aln.set_tag('RF', 'Y', value_type='Z')
                            continue
                       
                        any_change = True
                        # check
                        for n, (aln, flg) in enumerate(alns):
                            if (flg in primary_flags or aln.is_supplementary):
                                ps_aln = aln
                                if aln.get_tag('AS') == best_aln.get_tag('AS'):
                                    assert n == 0
                        dont_change = False
                        if len(alns) >= 2:
                            second_aln = alns[1][0]
                            if second_aln.get_tag('AS') == best_aln.get_tag('AS'):
                                dont_change = True 
                        
                        if not dont_change:
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
                            if designate_mapq != -1 and best_aln.mapping_quality < designate_mapq:
                                best_aln.mapping_quality = designate_mapq
                                best_aln.set_tag('RF', 'Y', value_type='Z')
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
                                output_alns.append(aln)
                    else:
                        for alns in output_aln_list:
                            for aln in alns:
                                output_alns.append(aln)
                else:
                    cached_aln.sort(key=lambda x: x[1])
                    for aln, flg in cached_aln:
                        if recalculate_all:
                            recalculate_score(aln, flg, ont_pos_index, ont_methylation_array)
                        output_alns.append(aln)
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
                    ont_pos_index, ont_methylation_array = get_5mc_sites_from_read(read_seq, mm_tag, ml_tag, prob_cutoff)
                    recalculate_score(aln, flg, ont_pos_index, ont_methylation_array)
                else:
                    pass
                    # print('{}\t{}\t{}\t{}\t{}\t{}\t{}\twhatever'.format(
                    #     last_read, aln.reference_name, aln.query_alignment_length, flg, aln.mapping_quality, aln.pos, aln.get_tag('AS')))
                output_alns.append(aln)
        
        return output_alns

    base_bam = os.path.basename(bam)
    out_bam = '{}.penalty{}.{}'.format(os.path.splitext(base_bam)[0], 
                                       penalty, 'methy_filtered.bam')
    
    format_options = [b'filter=!flag.unmap']
    def iter_read_groups(fin):
        last_name = None
        buf = []
        for aln in fin:
            name = aln.query_name
            if last_name is None:
                last_name = name
            if name != last_name:
                yield last_name, buf, aln.flag
                buf = []
                last_name = name
            buf.append((aln, aln.flag))
        if buf:
            yield last_name, buf, aln.flag

    def worker(group):
        read_name, alns, flag = group
        return analyze_cached_aln( alns, flag)


    with pysam.AlignmentFile(bam, format_options=format_options,
                             threads=threads) as fin, \
         pysam.AlignmentFile(out_bam, 'wb', template=fin,
                             threads=threads) as fout: 

         for group in iter_read_groups(fin):
            result = worker(group)
            for aln in result:
                fout.write(aln)
      
    return out_bam

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', help='input fasta file')
    parser.add_argument('bed', help='bed file output by pb-CpG-tools')
    parser.add_argument('bam', help='input bam file of ONT read mapping results')
    parser.add_argument('--penalty', type=int, default=2, help='penalty for inconsistent 5mC, default: %(default)s')
    parser.add_argument('--prob_cutoff', type=int, default=128, help='probability cutoff, default: %(default)s')
    parser.add_argument('--designate_mapq', type=int, default=60, help='designate MAPQ for the best alignments. Set the value to -1 to disable the function, default: %(default)s')
    parser.add_argument('--', default=False, action='store_true', help='recalculate scores for all alignments although there may not be any secondary alignments, default: %(default)s')
    parser.add_argument('--threads', type=int, default=8, help='number of threads for bam reading, default: %(default)s')
    args = parser.parse_args()
    
    if not -1 <= args.designate_mapq <= 60:
        raise Exception('--designate_mapq should be within the range [-1, 60]')
    
    fa_dict = parse_fasta(args.fasta)
    hifi_methylation_dict = parse_bedgraph(args.bed)
    parse_bam(args.bam, fa_dict, hifi_methylation_dict, args.penalty, args.prob_cutoff, args.designate_mapq, args.recalculate_all, args.threads)


if __name__ == '__main__':
    main()

