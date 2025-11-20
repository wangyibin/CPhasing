#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
parse gfa 
"""

import argparse
import logging
import os
import os.path as op
import sys

import re

import networkx as nx
import pandas as pd
import numpy as np

from collections import defaultdict, deque  
from pathlib import Path 
from pyranges import PyRanges
from concurrent.futures import ThreadPoolExecutor, as_completed

from .utilities import xopen 

logger = logging.getLogger(__name__)

_overlap_re = re.compile(r'^(\d+)M$')  

class Gfa:
    def __init__(self, gfa):
        self.file = gfa
        self.filename = Path(gfa).name
        self.segments = {}      # name -> length
        self.links = []         # list of dict(from, to, from_orient, to_orient, overlap_len, raw_overlap)
        self._graph = None
        self.parse()

    def parse(self):
        segments = {}
        segment_seqs = {}
        links = []
        with xopen(self.file) as f:
            for line in f:
                if not line.strip() or line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                t = parts[0]
                if t == 'S':
                    # GFA1: S <name> <sequence> [opt TAGs]
                    name = parts[1]
                    seq = parts[2]
                    length = None
           
                    for tag in parts[3:]:
                        if tag.startswith('LN:i:'):
                            try:
                                length = int(tag.split(':')[-1])
                            except:
                                pass
                            break
                    if length is None:
                        length = 0 if seq == '*' else len(seq)
                    segments[name] = length
                    if seq != '*':
                        segment_seqs[name] = seq
                elif t == 'L':
                    # L <from> <from_orient> <to> <to_orient> <overlap>
                    if len(parts) < 6:
                        continue
                    frm, fo, to, to_o, ov = parts[1:6]
                    ov_len = self._parse_overlap(ov)
                    links.append({
                        'from': frm,
                        'from_orient': fo,
                        'to': to,
                        'to_orient': to_o,
                        'overlap': ov,
                        'ov_len': ov_len
                    })
                else:
                    continue
        self.segments = segments
        self.seqment_seqs = segment_seqs
        if len(self.seqment_seqs) == 0:
            logger.warning("There is no segment sequences in the GFA file.")
        
        self.links = links
        self._graph = None

    @staticmethod
    def _parse_overlap(ov):
       
        m = _overlap_re.match(ov)
        if m:
            return int(m.group(1))
        num = ''
        for ch in ov:
            if ch.isdigit():
                num += ch
            else:
                break
        return int(num) if num else 0

    def to_graph(self):
        if self._graph is not None:
            return self._graph
        G = nx.DiGraph()
        for s, ln in self.segments.items():
            G.add_node(s, length=ln)
        for e in self.links:
            if e['from'] in G and e['to'] in G:
                G.add_edge(e['from'], e['to'],
                           ov_len=e['ov_len'],
                           overlap=e['overlap'],
                           from_orient=e['from_orient'],
                           to_orient=e['to_orient'])
        self._graph = G
        return G

    def remove_link_by_min_overlap(self, min_overlap=0):
        if min_overlap <= 0:
            return
        self.links = [e for e in self.links if e['ov_len'] >= min_overlap]
        self._graph = None

   
    def find_bubbles(self,
                     max_depth=50,
                     max_nodes_per_branch=200,
                     max_path_bp=500000,
                     require_all_branches=False):
        
        G = self.to_graph()
        bubbles = {}
        for s in G.nodes():
            succs = list(G.successors(s))
            if len(succs) < 2:
                continue
       
            sink_paths = defaultdict(list)
            for i in range(len(succs)):
                for j in range(i + 1, len(succs)):
                    u, v = succs[i], succs[j]
                    meet_node, path_u, path_v = self._pair_branch_meet(
                        G, s, u, v,
                        max_depth=max_depth,
                        max_nodes=max_nodes_per_branch,
                        max_bp=max_path_bp
                    )
                    if meet_node:
                        sink_paths[meet_node].append(path_u)
                        sink_paths[meet_node].append(path_v)
            if not sink_paths:
                continue
           
            for sink, paths in sink_paths.items():
                uniq = []
                seen_tuple = set()
                for p in paths:
                    t = tuple(p)
                    if t not in seen_tuple:
                        seen_tuple.add(t)
                        uniq.append(p)
                if require_all_branches and len(uniq) < len(succs):
                    continue
                key = (s, sink)
                if key not in bubbles:
                    bubbles[key] = {
                        'source': s,
                        'sink': sink,
                        'paths': uniq
                    }
                else:
                    for p in uniq:
                        t = tuple(p)
                        if all(tuple(x) != t for x in bubbles[key]['paths']):
                            bubbles[key]['paths'].append(p)


        result = []
        for (s, t), info in bubbles.items():
            second_nodes = set()
            for p in info['paths']:
                if len(p) >= 2:
                    second_nodes.add(p[1])
            info['n_branches'] = len(second_nodes)
            result.append(info)
        return result

    def _pair_branch_meet(self, G, source, u, v,
                          max_depth=50,
                          max_nodes=200,
                          max_bp=500000):
       
        def limited_bfs(start):
            paths = {start: [source, start]}
            q = deque([(start, 0, self.segments.get(start, 0))])
            while q:
                node, depth, acc_bp = q.popleft()
                if depth >= max_depth:
                    continue
                for nxt in G.successors(node):
                    if nxt == source:
                        continue  
                    edge_len = self.segments.get(nxt, 0)
                    new_bp = acc_bp + edge_len
                    if new_bp > max_bp:
                        continue
                    if nxt not in paths:
                        prev_path = paths[node][:]
                        paths[nxt] = prev_path + [nxt]
                        if len(paths) > max_nodes:
                            return paths 
                        q.append((nxt, depth + 1, new_bp))
            return paths

        paths_u = limited_bfs(u)
        paths_v = limited_bfs(v)

        candidates = set(paths_u.keys()) & set(paths_v.keys())
        candidates.discard(u)
        candidates.discard(v)
        candidates.discard(source)
        if not candidates:
            return None, None, None

        best = None
        best_len_sum = 10**12
        for c in candidates:
            lu = len(paths_u[c])
            lv = len(paths_v[c])
            if lu + lv < best_len_sum:
                best = c
                best_len_sum = lu + lv
        return best, paths_u[best], paths_v[best]


    def bubbles_to_dataframe(self, bubbles):
        rows = []
        for b in bubbles:
            for p in b['paths']:
                rows.append({
                    'source': b['source'],
                    'sink': b['sink'],
                    'n_branches': b['n_branches'],
                    'path': ','.join(p),
                    'path_len_nodes': len(p),
                    'path_len_bp': sum(self.segments.get(x, 0) for x in p)
                })
        return pd.DataFrame(rows)

    def filter_bubbles_two_utg(self,
                               bubbles,
                               utg_prefix='utg',
                               count_internal=True,
                               require_two_branches=True):
        res = []
        for b in bubbles:
            if require_two_branches and b.get('n_branches', None) != 2:
                continue
            all_nodes = set()
            internal_nodes = set()
            for p in b['paths']:
                all_nodes.update(p)
                if len(p) > 2:
                    internal_nodes.update(p[1:-1])  
            target_nodes = internal_nodes if count_internal else all_nodes
            utg_nodes = [n for n in target_nodes if n.startswith(utg_prefix)]
            if len(utg_nodes) == 2:
                b_copy = dict(b)
                b_copy['utg_nodes'] = utg_nodes
                res.append(b_copy)
        return res

    def bubbles_two_utg_to_dataframe(self, bubbles, **kwargs):

        filt = self.filter_bubbles_two_utg(bubbles, **kwargs)
        rows = []
        for b in filt:
            rows.append({
                'source': b['source'],
                'sink': b['sink'],
                'n_branches': b['n_branches'],
                'utg_nodes': ';'.join(b.get('utg_nodes', [])),
                'n_paths': len(b['paths'])
            })
        return pd.DataFrame(rows)

    def filter_bubbles_two_internal_utg_pattern(self,
                                                bubbles,
                                                utg_prefix='utg',
                                                exact_path_len=3,
                                                require_degree=True):

   
        G = self.to_graph()
        res = []
        for b in bubbles:
            if b.get('n_branches', 0) != 2:
                continue
            paths = b['paths']
            if len(paths) != 2:
                continue
            p1, p2 = paths
            if exact_path_len and (len(p1) != exact_path_len or len(p2) != exact_path_len):
                continue
            if not (p1[0] == p2[0] and p1[-1] == p2[-1]):
                continue
            source = p1[0]
            sink = p1[-1]
            internal1 = p1[1:-1]
            internal2 = p2[1:-1]
        
            if len(internal1) != 1 or len(internal2) != 1:
                continue
            n1, n2 = internal1[0], internal2[0]
            if n1 == n2:
                continue

            if not (n1.startswith(utg_prefix) and n2.startswith(utg_prefix)):
                continue
            if require_degree:
              
                preds1 = set(G.predecessors(n1)) - {source}
                preds2 = set(G.predecessors(n2)) - {source}
        
                succs1 = set(G.successors(n1)) - {sink}
                succs2 = set(G.successors(n2)) - {sink}
                if preds1 or preds2 or succs1 or succs2:
                    continue
            b_copy = dict(b)
            b_copy['internal_nodes'] = [n1, n2]
            res.append(b_copy)
        return res

    def bubbles_two_internal_utg_to_dataframe(self, bubbles, **kwargs):
   
        filt = self.filter_bubbles_two_internal_utg_pattern(bubbles, **kwargs)
        rows = []
        for b in filt:
            rows.append({
                'source': b['source'],
                'sink': b['sink'],
                'internal_utg1': b['internal_nodes'][0],
                'internal_utg2': b['internal_nodes'][1],
                'n_paths': len(b['paths'])
            })
        return pd.DataFrame(rows)

    def overlap_to_bed(self, invert=False):
        res = []
        for e in self.links:
            contig1 = e['from']
            strand1 = e['from_orient']
            contig2 = e['to']
            strand2 = e['to_orient']
            overlap = e['ov_len']
            length1 = self.segments[contig1]
            length2 = self.segments[contig2]
            if strand1 == "+":
                res.append((contig1, length1 - overlap, length1))
            else:
                res.append((contig1, 0, overlap))
        
            if strand2 == "-":
                res.append((contig2, length2 - overlap, length2))
            else:
                res.append((contig2, 0, overlap))

        df = pd.DataFrame(res).sort_values([0, 1])
        df.columns = ['Chromosome', 'Start', 'End']
        df = PyRanges(df).merge().df

        if invert:
            chrom_df = pd.DataFrame({
                'Chromosome': list(self.segments.keys()),
                'Start': [0] * len(self.segments),
                'End': list(self.segments.values())
            })
            chrom_df = PyRanges(chrom_df)
            bed_df = PyRanges(df)
            inv = chrom_df.subtract(bed_df).df
            return inv

        return df
    


    @staticmethod
    def _revcomp(seq: str) -> str:
        comp = str.maketrans('ACGTRYMKBDHVNacgtrymkbdhvn',
                             'TGCAYRKMVHDBNtgcayrkmvhdbn')
        return seq.translate(comp)[::-1]

    def _find_link_between(self, a, b):
     
        for e in self.links:
            if e['from'] == a and e['to'] == b:
                return e, 1
            if e['from'] == b and e['to'] == a:
                return e, -1
        return None, 0

    def _internal_unique_sequence(self, source, internal, sink):
        
        if internal not in self.seqment_seqs:
            raise ValueError(f"The sequence missing {internal} (S line might be '*')")

        seq = self.seqment_seqs[internal]
        seq_oriented = seq
        oriented = False 

        left_ov = 0
        e1, dir1 = self._find_link_between(source, internal)
        if e1:
            if dir1 == 1:
  
                left_ov = e1['ov_len']
                internal_orient = e1['to_orient']
            else:
     
                left_ov = e1['ov_len']
                internal_orient = e1['from_orient']
            if internal_orient == '-' and not oriented:
                seq_oriented = self._revcomp(seq_oriented)
                oriented = True


        right_ov = 0
        e2, dir2 = self._find_link_between(internal, sink)
        if e2:
            if dir2 == 1:
  
                internal_orient2 = e2['from_orient']
                right_ov = e2['ov_len']
            else:

                internal_orient2 = e2['to_orient']
                right_ov = e2['ov_len']
         
            if internal_orient2 == '-' and not oriented:
                seq_oriented = self._revcomp(seq_oriented)
                oriented = True

        left_cut = left_ov
        right_cut = right_ov
        if left_cut + right_cut >= len(seq_oriented):

            left_cut = min(left_cut, max(0, len(seq_oriented) - 1))
            right_cut = 0
        unique_seq = seq_oriented[left_cut: len(seq_oriented) - right_cut if right_cut > 0 else len(seq_oriented)]
        return unique_seq, {
            'left_overlap': left_ov,
            'right_overlap': right_ov,
            'final_length': len(unique_seq),
            'oriented': oriented
        }

    def extract_bubble_internal_unique_sequences(self, bubble):

        if 'internal_nodes' in bubble:
            internals = bubble['internal_nodes']
        else:
 
            paths = bubble['paths']
            mids = set()
            for p in paths:
                if len(p) >= 3:
                    mids.add(p[1])
            if len(mids) != 2:
                raise ValueError("Two internal nodes cannot be automatically determined (bubble['internal_nodes'] is required).")
            internals = list(mids)

        source = bubble['source']
        sink = bubble['sink']
        result = {}
        for internal in internals:
            seq_trim, meta = self._internal_unique_sequence(source, internal, sink)
            result[internal] = {'seq': seq_trim, 'meta': meta}
        return result

    def write_bubble_internal_unique_fasta(self, bubble, out_fa):
        data = self.extract_bubble_internal_unique_sequences(bubble)
        with open(out_fa, 'w') as fh:
            for name, info in data.items():
                seq = info['seq']
                meta = info['meta']
                desc = f"leftOv={meta['left_overlap']} rightOv={meta['right_overlap']} len={meta['final_length']}"
                fh.write(f">{name} {desc}\n")
                
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i+80] + '\n')
        return out_fa
    
    def _revcomp2(self, s: str) -> str:
        tbl = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
        return s.translate(tbl)[::-1]

    def _parse_edlib_cigar(self, cigar: str):
     
        m = x = ins = dele = 0
        num = []
        for ch in cigar or "":
            if ch.isdigit():
                num.append(ch)
                continue
            n = int(''.join(num)) if num else 1
            num = []
            if ch == '=':
                m += n
            elif ch == 'X':
                x += n
            elif ch == 'I':
                ins += n
            elif ch == 'D':
                dele += n
            elif ch == 'M': 
                m += n
        alen = m + x + ins + dele
        q_consume = m + x + ins
        t_consume = m + x + dele
        return m, x, ins, dele, alen, q_consume, t_consume
    
    def _cigar_qt_consume(self, cigar: str):
        q_consume = t_consume = alen = 0
        num = []
        for ch in cigar or "":
            if ch.isdigit():
                num.append(ch)
                continue
            n = int(''.join(num)) if num else 1
            num = []
            if ch in ('M', '=', 'X'):
                q_consume += n
                t_consume += n
                alen += n  
            elif ch == 'I':
                q_consume += n
                alen += n
            elif ch == 'D':
                t_consume += n
                alen += n
        return q_consume, t_consume, alen
    
    def align_internal_utgs(self, bubble, output_prefix="bubble_utg_align",
                            method="wfmash", minimap2="minimap2",
                            biopython_mode="globalxx", mm2_threads=1,
                            ps_match=2, ps_mismatch=-3, ps_open=5, ps_extend=2):
     
        fasta = f"{output_prefix}.fa"
        self.write_bubble_internal_unique_fasta(bubble, fasta)
        internals = list(self.extract_bubble_internal_unique_sequences(bubble).keys())
        if len(internals) != 2:
            raise ValueError("Need two internal utg")
        a, b2 = internals
        if len(internals) != 2:
            raise ValueError("Need two internal utg")
        if method == "minimap2":
            paf = f"{output_prefix}.paf"

            cmd = f"{minimap2} -t {mm2_threads} -c -k 19 -w 19 -DP {fasta} {fasta} > {paf}"
            import subprocess
            subprocess.run(cmd, shell=True, check=True)
            return {'fasta': fasta, 'result': paf, 'cmd': cmd}
        
        elif method == "wfmash":
            paf = f"{output_prefix}.paf"
            cmd = f"wfmash -t {mm2_threads} -n 1 -X -N {fasta} {fasta} > {paf}"
            import subprocess
            subprocess.run(cmd, shell=True, check=True)
            return {'fasta': fasta, 'result': paf, 'cmd': cmd}
        elif method == "edlib":
            try:
                import edlib
            except ImportError as e:
                raise ImportError("Please intall edlib: pip install edlib") from e

            data = self.extract_bubble_internal_unique_sequences(bubble)
            a, b2 = internals
            s1 = data[a]['seq']
            s2 = data[b2]['seq']
            qlen, tlen = len(s1), len(s2)

            def parse_ext_cigar(cigar: str):
                m = x = ins = dele = 0
                num = []
                for ch in cigar or "":
                    if ch.isdigit():
                        num.append(ch)
                        continue
                    n = int(''.join(num)) if num else 1
                    num = []
                    if ch == '=':
                        m += n
                    elif ch == 'X':
                        x += n
                    elif ch == 'I':
                        ins += n
                    elif ch == 'D':
                        dele += n
                    elif ch == 'M': 
                        m += n
                alen = m + x + ins + dele
                q_consume = m + x + ins
                t_consume = m + x + dele
                return m, x, ins, dele, alen, q_consume, t_consume

            def ext_to_cg(cigar: str):
                out = []
                num = []
                last_op = None
                run = 0
                def flush(op, n):
                    if n > 0:
                        out.append(f"{n}{op}")
                for ch in cigar or "":
                    if ch.isdigit():
                        num.append(ch)
                        continue
                    n = int(''.join(num)) if num else 1
                    num = []
                    op = 'M' if ch in ('=', 'X', 'M') else ch
                    if op == last_op:
                        run += n
                    else:
                        if last_op is not None:
                            flush(last_op, run)
                        last_op, run = op, n
                if last_op is not None:
                    flush(last_op, run)
                return ''.join(out)

            def run_one(target: str, strand: str):
                res = edlib.align(s1, target, mode="HW", task="path")
                cigar = res.get('cigar') or ""
                # edlib.locations: [(start, end)] 0-based, inclusive
                locs = res.get('locations') or []
                if locs and isinstance(locs[0], (list, tuple)) and len(locs[0]) == 2:
                    ts_incl, te_incl = int(locs[0][0]), int(locs[0][1])
                else:
                    ts_incl, te_incl = 0, len(target) - 1
                m, x, ins, dele, alen, q_consume, t_consume = parse_ext_cigar(cigar)
                nmatch = m
                ident = (nmatch / alen) if alen > 0 else 0.0

                qstart, qend = 0, q_consume
                if strand == '+':
                    tstart, tend = ts_incl, te_incl + 1
                else:
                    tstart = tlen - 1 - te_incl
                    tend   = (tlen - 1 - ts_incl) + 1
                return {
                    'identity': ident, 'nmatch': nmatch, 'alen': alen,
                    'qstart': qstart, 'qend': qend,
                    'tstart': tstart, 'tend': tend,
                    'strand': strand, 'cigar_ext': cigar
                }

            info_f = run_one(s2, '+')
            info_r = run_one(self._revcomp2(s2), '-')
            best = info_f if (info_f['identity'], info_f['alen']) >= (info_r['identity'], info_r['alen']) else info_r

            paf = f"{output_prefix}.paf"
            mapq = 255
            cg = ext_to_cg(best['cigar_ext']) if best['cigar_ext'] else ""
            with open(paf, 'w') as fh:
                fields = [
                    a, str(qlen), str(best['qstart']), str(best['qend']),
                    best['strand'],
                    b2, str(tlen), str(best['tstart']), str(best['tend']),
                    str(best['nmatch']), str(best['alen']), str(mapq)
                ]
                line = '\t'.join(fields)
                if cg:
                    line += f"\tcg:Z:{cg}\tNM:i:{best['alen'] - best['nmatch']}"
                fh.write(line + '\n')

            return {'fasta': fasta, 'result': paf, 'cmd': 'edlib:HW'}
        
        elif method == "biopython":
            try:
                from Bio import pairwise2
            except ImportError:
                raise ImportError("Can not import biopython for method=biopython")
            seqs = {}
            for line in open(fasta):
                if line.startswith('>'):
                    cur = line[1:].split()[0]
                    seqs[cur] = []
                else:
                    seqs[cur].append(line.strip())
            seqs = {k: ''.join(v) for k, v in seqs.items()}
            a, b = internals
            s1 = seqs[a]
            s2 = seqs[b]
            if biopython_mode == "globalxx":
                alns = pairwise2.align.globalxx(s1, s2, one_alignment_only=True)
            else:
                alns = pairwise2.align.localxx(s1, s2, one_alignment_only=True)
            out_txt = f"{output_prefix}.biopy.aln.txt"
            with open(out_txt, 'w') as fh:
                if alns:
                    aln = alns[0]
                    fh.write(f"{a} vs {b}\nScore={aln.score}\n{aln.seqA}\n{aln.seqB}\n")
                else:
                    fh.write("No alignment produced\n")
            return {'fasta': fasta, 'result': out_txt}
        else:
            raise ValueError("method only support for minimap2 or biopython")

    def _ensure_dir(self, d):
        if d and not os.path.isdir(d):
            os.makedirs(d, exist_ok=True)

    def _parse_pair_paf(self, paf_path, a, b):
  
        best = None
        if not os.path.isfile(paf_path):
            return best
        def _to_int(x, default=0):
            try:
                return int(x)
            except Exception:
                return default
        with open(paf_path) as fh:
            for line in fh:
                if not line.strip() or line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 12:
                    continue
                qn, qlen, qs, qe, strand, tn, tlen, ts, te, nmatch, alen, mapq = parts[:12]
                if qn == tn:
                    continue
                if not ({qn, tn} == {a, b}):
                    continue
                qlen = _to_int(qlen)
                qs, qe = _to_int(qs), _to_int(qe)
                tlen = _to_int(tlen)
                ts, te = _to_int(ts), _to_int(te)
                nmatch = _to_int(nmatch)
                alen = _to_int(alen, default=max(qe-qs, te-ts))
                mapq = _to_int(mapq, default=0)
                ident = (nmatch / alen) if alen > 0 else 0.0
                cg = None
                for tag in parts[12:]:
                    if tag.startswith('cg:Z:'):
                        cg = tag[5:]
                        break
                rec = {
                    'qname': qn, 'tname': tn, 'strand': strand,
                    'qlen': qlen, 'qstart': qs, 'qend': qe,
                    'tlen': tlen, 'tstart': ts, 'tend': te,
                    'nmatch': nmatch, 'alen': alen, 'mapq': mapq,
                    'identity': ident, 'cigar': cg
                }
                if (best is None) or (rec['identity'] > best['identity']) or \
                   (rec['identity'] == best['identity'] and rec['alen'] > best['alen']):
                    best = rec
        return best

    def _parse_edlib_cigar(self, cigar: str):
        m = x = ins = dele = 0
        num = []
        for ch in cigar:
            if ch.isdigit():
                num.append(ch)
                continue
            n = int(''.join(num)) if num else 1
            num = []
            if ch == '=':
                m += n
            elif ch == 'X':
                x += n
            elif ch == 'I':
                ins += n
            elif ch == 'D':
                dele += n
            elif ch == 'M':  
                m += n
        alen = m + x + ins + dele
        return m, x, ins, dele, alen


    def align_bubbles_and_collect(self, bubbles,
                                  out_dir="bubble_align_all",
                                  utg_prefix='utg',
                                  exact_path_len=3,
                                  require_degree=False,
                                  method="wfmash",
                                  minimap2="minimap2",
                                  biopython_mode="globalxx",
                                  tsv_path=None,
                                  n_jobs=20,             
                                  mm2_threads=1):       

        targets = self.filter_bubbles_two_internal_utg_pattern(
            bubbles, utg_prefix=utg_prefix,
            exact_path_len=exact_path_len,
            require_degree=require_degree
        )
        self._ensure_dir(out_dir)

        def do_one(idx, b):
            bid = f"b{idx}"
            internals = b['internal_nodes']
            data = self.extract_bubble_internal_unique_sequences(b)
            a, b2 = internals[0], internals[1]
            len_a = len(data[a]['seq'])
            len_b = len(data[b2]['seq'])
            meta_a = data[a]['meta']
            meta_b = data[b2]['meta']
            prefix = os.path.join(out_dir, f"{bid}_{b['source']}_{b['sink']}")
            try:
                res = self.align_internal_utgs(
                    b, output_prefix=prefix, method=method,
                    minimap2=minimap2, biopython_mode=biopython_mode,
                    mm2_threads=mm2_threads
                )
                fasta_path = res['fasta']
                out_path = res['result']
                if method == "minimap2":
                    paf_best = self._parse_pair_paf(out_path, a, b2)
                    if paf_best is None:
                        row = {
                            'bubble_id': bid, 'source': b['source'], 'sink': b['sink'],
                            'int1': a, 'int2': b2, 'len1': len_a, 'len2': len_b,
                            'leftOv1': meta_a['left_overlap'], 'rightOv1': meta_a['right_overlap'],
                            'leftOv2': meta_b['left_overlap'], 'rightOv2': meta_b['right_overlap'],
                            'method': method, 'fasta': fasta_path, 'paf': out_path,
                            'identity': np.nan, 'alen': 0, 'nmatch': 0, 'mapq': np.nan,
                            'strand': None, 'qstart': None, 'qend': None, 'tstart': None, 'tend': None, 'cigar': None
                        }
                    else:
                        row = {
                            'bubble_id': bid, 'source': b['source'], 'sink': b['sink'],
                            'int1': a, 'int2': b2, 'len1': len_a, 'len2': len_b,
                            'leftOv1': meta_a['left_overlap'], 'rightOv1': meta_a['right_overlap'],
                            'leftOv2': meta_b['left_overlap'], 'rightOv2': meta_b['right_overlap'],
                            'method': method, 'fasta': fasta_path, 'paf': out_path,
                            'identity': paf_best['identity'], 'alen': paf_best['alen'],
                            'nmatch': paf_best['nmatch'], 'mapq': paf_best['mapq'],
                            'strand': paf_best['strand'],
                            'qstart': paf_best['qstart'], 'qend': paf_best['qend'],
                            'tstart': paf_best['tstart'], 'tend': paf_best['tend'],
                            'cigar': paf_best['cigar']
                        }
                elif method in ("edlib", "parasail"):
                    st = res.get('stats', {})
                    row = {
                        'bubble_id': bid, 'source': b['source'], 'sink': b['sink'],
                        'int1': a, 'int2': b2, 'len1': len_a, 'len2': len_b,
                        'leftOv1': meta_a['left_overlap'], 'rightOv1': meta_a['right_overlap'],
                        'leftOv2': meta_b['left_overlap'], 'rightOv2': meta_b['right_overlap'],
                        'method': method, 'fasta': fasta_path, 'paf_or_txt': out_path,
                        'identity': st.get('identity', np.nan), 'alen': st.get('alen', 0),
                        'nmatch': st.get('nmatch', 0), 'mapq': np.nan,
                        'strand': st.get('strand', '+'),
                        'qstart': st.get('qstart'), 'qend': st.get('qend'),
                        'tstart': st.get('tstart'), 'tend': st.get('tend'),
                        'cigar': st.get('cigar')
                    }
                else:
                    row = {
                        'bubble_id': bid, 'source': b['source'], 'sink': b['sink'],
                        'int1': a, 'int2': b2, 'len1': len_a, 'len2': len_b,
                        'leftOv1': meta_a['left_overlap'], 'rightOv1': meta_a['right_overlap'],
                        'leftOv2': meta_b['left_overlap'], 'rightOv2': meta_b['right_overlap'],
                        'method': method, 'fasta': fasta_path, 'aln_txt': out_path
                    }
            except Exception as e:
                logger.warning(f"bubble {bid} failed align: {e}")
                row = {
                    'bubble_id': bid, 'source': b['source'], 'sink': b['sink'],
                    'int1': a, 'int2': b2, 'len1': len_a, 'len2': len_b,
                    'leftOv1': meta_a['left_overlap'], 'rightOv1': meta_a['right_overlap'],
                    'leftOv2': meta_b['left_overlap'], 'rightOv2': meta_b['right_overlap'],
                    'method': method, 'fasta': None, 'paf': None,
                    'identity': np.nan, 'alen': 0, 'nmatch': 0, 'mapq': np.nan,
                    'strand': None, 'qstart': None, 'qend': None, 'tstart': None, 'tend': None, 'cigar': None
                }
            return idx, row

        rows = [None] * len(targets)
        with ThreadPoolExecutor(max_workers=max(1, int(n_jobs))) as ex:
            futs = {ex.submit(do_one, idx, b): idx for idx, b in enumerate(targets)}
            for fut in as_completed(futs):
                idx, row = fut.result()
                rows[idx] = row

        df = pd.DataFrame([r for r in rows if r is not None])
        if tsv_path:
            df.to_csv(tsv_path, sep='\t', index=False)
        return df    

    def write(self, out_handle):
        for name, ln in self.segments.items():
            out_handle.write(f"S\t{name}\t*\tLN:i:{ln}\n")
        for e in self.links:
            length1 = self.segments[e['from']]
            length2 = self.segments[e['to']]
            out_handle.write(f"L\t{e['from']}\t{e['from_orient']}\t{e['to']}\t{e['to_orient']}\t{e['overlap']}\tL1:i:{length1}\tL2:i:{length2}\n")
