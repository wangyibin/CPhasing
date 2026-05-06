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

_overlap_re = re.compile(r"^(\d+)M$")


class Gfa:
    def __init__(self, gfa):
        self.file = gfa
        self.filename = Path(gfa).name
        self.segments = {}  # name -> length
        self.links = []  # list of dict(from, to, from_orient, to_orient, overlap_len, raw_overlap)
        self._graph = None
        self.parse()

    def parse(self):
        segments = {}
        segment_seqs = {}
        segment_tags = {}
        links = []
        with xopen(self.file) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                t = parts[0]
                if t == "S":
                    # GFA1: S <name> <sequence> [opt TAGs]
                    name = parts[1]
                    seq = parts[2]
                    length = None
                    tags = parts[3:]
                    for tag in tags:
                        if tag.startswith("LN:i:"):
                            try:
                                length = int(tag.split(":")[-1])
                            except:
                                pass
                            break
                    if length is None:
                        length = 0 if seq == "*" else len(seq)
                    segments[name] = length
                    segment_tags[name] = tags
                   
                    if seq != "*":
                        segment_seqs[name] = seq
                elif t == "L":
                    # L <from> <from_orient> <to> <to_orient> <overlap>
                    if len(parts) < 6:
                        continue
                    frm, fo, to, to_o, ov = parts[1:6]
                    ov_len = self._parse_overlap(ov)
                    links.append(
                        {
                            "from": frm,
                            "from_orient": fo,
                            "to": to,
                            "to_orient": to_o,
                            "overlap": ov,
                            "ov_len": ov_len,
                        }
                    )
                else:
                    continue
        self.segments = segments
        self.seqment_seqs = segment_seqs
        self.segment_tags = segment_tags
        if len(self.seqment_seqs) == 0:
            logger.warning("There is no segment sequences in the GFA file.")

        self.links = links
        self._graph = None

    @staticmethod
    def _parse_overlap(ov):
        m = _overlap_re.match(ov)
        if m:
            return int(m.group(1))
        num = ""
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
            if e["from"] in G and e["to"] in G:
                G.add_edge(
                    e["from"],
                    e["to"],
                    ov_len=e["ov_len"],
                    overlap=e["overlap"],
                    from_orient=e["from_orient"],
                    to_orient=e["to_orient"],
                )
        self._graph = G
        return G

    def remove_link_by_min_overlap(self, min_overlap=0):
        if min_overlap <= 0:
            return
        self.links = [e for e in self.links if e["ov_len"] >= min_overlap]
        self._graph = None

    def find_bubbles(
        self,
        max_depth=50,
        max_nodes_per_branch=200,
        max_path_bp=500000,
        require_all_branches=False,
    ):
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
                        G,
                        s,
                        u,
                        v,
                        max_depth=max_depth,
                        max_nodes=max_nodes_per_branch,
                        max_bp=max_path_bp,
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
                    bubbles[key] = {"source": s, "sink": sink, "paths": uniq}
                else:
                    for p in uniq:
                        t = tuple(p)
                        if all(tuple(x) != t for x in bubbles[key]["paths"]):
                            bubbles[key]["paths"].append(p)

        result = []
        for (s, t), info in bubbles.items():
            second_nodes = set()
            for p in info["paths"]:
                if len(p) >= 2:
                    second_nodes.add(p[1])
            info["n_branches"] = len(second_nodes)
            result.append(info)
        return result

    def _pair_branch_meet(
        self, G, source, u, v, max_depth=50, max_nodes=200, max_bp=500000
    ):
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
            for p in b["paths"]:
                rows.append(
                    {
                        "source": b["source"],
                        "sink": b["sink"],
                        "n_branches": b["n_branches"],
                        "path": ",".join(p),
                        "path_len_nodes": len(p),
                        "path_len_bp": sum(self.segments.get(x, 0) for x in p),
                    }
                )
        return pd.DataFrame(rows)

    def filter_bubbles_two_utg(
        self, bubbles, utg_prefix="utg", count_internal=True, require_two_branches=True
    ):
        res = []
        for b in bubbles:
            if require_two_branches and b.get("n_branches", None) != 2:
                continue
            all_nodes = set()
            internal_nodes = set()
            for p in b["paths"]:
                all_nodes.update(p)
                if len(p) > 2:
                    internal_nodes.update(p[1:-1])
            target_nodes = internal_nodes if count_internal else all_nodes
            utg_nodes = [n for n in target_nodes if n.startswith(utg_prefix)]
            if len(utg_nodes) == 2:
                b_copy = dict(b)
                b_copy["utg_nodes"] = utg_nodes
                res.append(b_copy)
        return res

    def bubbles_two_utg_to_dataframe(self, bubbles, **kwargs):
        filt = self.filter_bubbles_two_utg(bubbles, **kwargs)
        rows = []
        for b in filt:
            rows.append(
                {
                    "source": b["source"],
                    "sink": b["sink"],
                    "n_branches": b["n_branches"],
                    "utg_nodes": ";".join(b.get("utg_nodes", [])),
                    "n_paths": len(b["paths"]),
                }
            )
        return pd.DataFrame(rows)

    def filter_bubbles_two_internal_utg_pattern(
        self, bubbles, utg_prefix="utg", exact_path_len=3, require_degree=True
    ):
        G = self.to_graph()
        res = []
        for b in bubbles:
            if b.get("n_branches", 0) != 2:
                continue
            paths = b["paths"]
            if len(paths) != 2:
                continue
            p1, p2 = paths
            if exact_path_len and (
                len(p1) != exact_path_len or len(p2) != exact_path_len
            ):
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
            b_copy["internal_nodes"] = [n1, n2]
            res.append(b_copy)
        return res

    def bubbles_two_internal_utg_to_dataframe(self, bubbles, **kwargs):
        filt = self.filter_bubbles_two_internal_utg_pattern(bubbles, **kwargs)
        rows = []
        for b in filt:
            rows.append(
                {
                    "source": b["source"],
                    "sink": b["sink"],
                    "internal_utg1": b["internal_nodes"][0],
                    "internal_utg2": b["internal_nodes"][1],
                    "n_paths": len(b["paths"]),
                }
            )
        return pd.DataFrame(rows)

    def overlap_to_bed(self, invert=False):
        chroms = []
        starts = []
        ends = []
        
        segments = self.segments
        for e in self.links:
            c1, s1 = e["from"], e["from_orient"]
            c2, s2 = e["to"], e["to_orient"]
            ov = e["ov_len"]
            
            l1 = segments[c1]
            l2 = segments[c2]
            
            chroms.extend((c1, c2))
            
            if s1 == "+":
                starts.append(l1 - ov)
                ends.append(l1)
            else:
                starts.append(0)
                ends.append(ov)

            if s2 == "-":
                starts.append(l2 - ov)
                ends.append(l2)
            else:
                starts.append(0)
                ends.append(ov)
        df = pd.DataFrame({"Chromosome": chroms, "Start": starts, "End": ends})
        df = PyRanges(df).merge().df

        if invert:
            chrom_df = pd.DataFrame(
                {
                    "Chromosome": list(segments.keys()),
                    "Start": 0,
                    "End": list(segments.values()),
                }
            )
            inv = PyRanges(chrom_df).subtract(PyRanges(df)).df
            return inv

        return df

    def calculate_overlap_percentage(self):
        overlap_df = self.overlap_to_bed()

        if overlap_df.empty:
            result_df = pd.DataFrame(
                list(self.segments.items()), columns=["Chromosome", "Length"]
            )
            result_df["Overlap_Length"] = 0
            result_df["Overlap_Percentage"] = 0.0
            return result_df

        overlap_df["Overlap_Length"] = overlap_df["End"] - overlap_df["Start"]

        overlap_sum = (
            overlap_df.groupby("Chromosome")["Overlap_Length"].sum().reset_index()
        )

        contig_lengths = pd.DataFrame(
            list(self.segments.items()), columns=["Chromosome", "Length"]
        )

        result_df = pd.merge(contig_lengths, overlap_sum, on="Chromosome", how="left")

        result_df["Overlap_Length"] = result_df["Overlap_Length"].fillna(0).astype(int)

        result_df["Overlap_Percentage"] = (
            result_df["Overlap_Length"] / result_df["Length"]
        ) * 100
        result_df["Overlap_Percentage"] = result_df["Overlap_Percentage"].clip(
            upper=100.0
        )

        return result_df

    def summary(self):
        lengths = list(self.segments.values())
        num_seqs = len(lengths)
        sum_len = sum(lengths) if lengths else 0
        min_len = min(lengths) if lengths else 0
        max_len = max(lengths) if lengths else 0
        avg_len = sum_len / num_seqs if num_seqs > 0 else 0

        lengths.sort(reverse=True)
        n50, L50 = 0, 0
        n90, L90 = 0, 0
        cumsum = 0
        for i, l in enumerate(lengths):
            cumsum += l
            if n50 == 0 and cumsum >= sum_len / 2:
                n50 = l
                L50 = i + 1
            if n90 == 0 and cumsum >= sum_len * 0.9:
                n90 = l
                L90 = i + 1

        link_num = len(self.links)

        ov_lens = [e["ov_len"] for e in self.links if e["ov_len"] > 0]
        overlap_seq_len = sum(ov_lens) if ov_lens else 0
        overlap_seq_median = np.median(ov_lens) if ov_lens else 0
        overlap_seq_avg = np.mean(ov_lens) if ov_lens else 0

        overlap_df = self.calculate_overlap_percentage()
        if not overlap_df.empty:
            cov_100 = (overlap_df["Overlap_Percentage"] >= 99.99).sum()
            cov_100_len = (
                overlap_df[overlap_df["Overlap_Percentage"] >= 99.99]["Length"]
            ).sum()
            cov_90 = (overlap_df["Overlap_Percentage"] >= 90).sum()
            cov_90_len = (
                overlap_df[overlap_df["Overlap_Percentage"] >= 90]["Length"]
            ).sum()
            cov_75 = (overlap_df["Overlap_Percentage"] >= 75).sum()
            cov_75_len = (
                overlap_df[overlap_df["Overlap_Percentage"] >= 75]["Length"]
            ).sum()
            cov_50 = (overlap_df["Overlap_Percentage"] >= 50).sum()
            cov_50_len = (
                overlap_df[overlap_df["Overlap_Percentage"] >= 50]["Length"]
            ).sum()
        else:
            cov_100 = cov_90 = cov_75 = cov_50 = 0

        try:
            G = self.to_graph()
            isolated_nodes = nx.number_of_isolates(G)
            dead_ends = sum(
                1 for n in G.nodes() if G.in_degree(n) == 0 or G.out_degree(n) == 0
            )
            n_components = (
                nx.number_weakly_connected_components(G) if num_seqs > 0 else 0
            )
        except Exception as e:
            logger.warning(f"Graph topology calculation failed: {e}")
            isolated_nodes = dead_ends = n_components = 0

        avg_degree = (link_num * 2 / num_seqs) if num_seqs > 0 else 0.0

        try:
            bubbles = self.find_bubbles()
            total_bubbles = len(bubbles)
            simple_bubbles = sum(1 for b in bubbles if b.get("n_branches", 0) == 2)
        except Exception as e:
            logger.warning(f"Bubble calculation failed: {e}")
            total_bubbles = simple_bubbles = 0

        summary_dict = {
            "num_seqs": num_seqs,
            "sum_len": sum_len,
            "min_len": min_len,
            "avg_len": round(avg_len, 2),
            "max_len": max_len,
            "N50": n50,
            "L50": L50,
            "N90": n90,
            "L90": L90,
            "link_num": link_num,
            "avg_degree": round(avg_degree, 2),
            "isolated_nodes": isolated_nodes,
            "dead_ends": dead_ends,
            "connected_components": n_components,
            "total_bubbles": total_bubbles,
            "simple_bubbles": simple_bubbles,
            "overlap_seq_len": overlap_seq_len,
            "overlap_seq_median": round(float(overlap_seq_median), 2),
            "overlap_seq_avg": round(float(overlap_seq_avg), 2),
            "overlap_coverage_ge100": int(cov_100),
            "overlap_coverage_ge100_len": int(cov_100_len),
            "overlap_coverage_ge90": int(cov_90),
            "overlap_coverage_ge90_len": int(cov_90_len),
            "overlap_coverage_ge75": int(cov_75),
            "overlap_coverage_ge75_len": int(cov_75_len),
            "overlap_coverage_ge50": int(cov_50),
            "overlap_coverage_ge50_len": int(cov_50_len),
        }

        return pd.Series(summary_dict)

    @staticmethod
    def _revcomp(seq: str) -> str:
        comp = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")
        return seq.translate(comp)[::-1]

    def _find_link_between(self, a, b):
        for e in self.links:
            if e["from"] == a and e["to"] == b:
                return e, 1
            if e["from"] == b and e["to"] == a:
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
                left_ov = e1["ov_len"]
                internal_orient = e1["to_orient"]
            else:
                left_ov = e1["ov_len"]
                internal_orient = e1["from_orient"]
            if internal_orient == "-" and not oriented:
                seq_oriented = self._revcomp(seq_oriented)
                oriented = True

        right_ov = 0
        e2, dir2 = self._find_link_between(internal, sink)
        if e2:
            if dir2 == 1:
                internal_orient2 = e2["from_orient"]
                right_ov = e2["ov_len"]
            else:
                internal_orient2 = e2["to_orient"]
                right_ov = e2["ov_len"]

            if internal_orient2 == "-" and not oriented:
                seq_oriented = self._revcomp(seq_oriented)
                oriented = True

        left_cut = left_ov
        right_cut = right_ov
        if left_cut + right_cut >= len(seq_oriented):
            left_cut = min(left_cut, max(0, len(seq_oriented) - 1))
            right_cut = 0
        unique_seq = seq_oriented[
            left_cut : len(seq_oriented) - right_cut
            if right_cut > 0
            else len(seq_oriented)
        ]
        return unique_seq, {
            "left_overlap": left_ov,
            "right_overlap": right_ov,
            "final_length": len(unique_seq),
            "oriented": oriented,
        }

    def extract_bubble_internal_unique_sequences(self, bubble):
        if "internal_nodes" in bubble:
            internals = bubble["internal_nodes"]
        else:
            paths = bubble["paths"]
            mids = set()
            for p in paths:
                if len(p) >= 3:
                    mids.add(p[1])
            if len(mids) != 2:
                raise ValueError(
                    "Two internal nodes cannot be automatically determined (bubble['internal_nodes'] is required)."
                )
            internals = list(mids)

        source = bubble["source"]
        sink = bubble["sink"]
        result = {}
        for internal in internals:
            seq_trim, meta = self._internal_unique_sequence(source, internal, sink)
            result[internal] = {"seq": seq_trim, "meta": meta}
        return result

    def write_bubble_internal_unique_fasta(self, bubble, out_fa):
        data = self.extract_bubble_internal_unique_sequences(bubble)
        with open(out_fa, "w") as fh:
            for name, info in data.items():
                seq = info["seq"]
                meta = info["meta"]
                desc = f"leftOv={meta['left_overlap']} rightOv={meta['right_overlap']} len={meta['final_length']}"
                fh.write(f">{name} {desc}\n")

                for i in range(0, len(seq), 80):
                    fh.write(seq[i : i + 80] + "\n")
        return out_fa

    def _revcomp2(self, s: str) -> str:
        tbl = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        return s.translate(tbl)[::-1]

    def _parse_edlib_cigar(self, cigar: str):
        m = x = ins = dele = 0
        num = []
        for ch in cigar or "":
            if ch.isdigit():
                num.append(ch)
                continue
            n = int("".join(num)) if num else 1
            num = []
            if ch == "=":
                m += n
            elif ch == "X":
                x += n
            elif ch == "I":
                ins += n
            elif ch == "D":
                dele += n
            elif ch == "M":
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
            n = int("".join(num)) if num else 1
            num = []
            if ch in ("M", "=", "X"):
                q_consume += n
                t_consume += n
                alen += n
            elif ch == "I":
                q_consume += n
                alen += n
            elif ch == "D":
                t_consume += n
                alen += n
        return q_consume, t_consume, alen

    def align_internal_utgs(
        self,
        bubble,
        output_prefix="bubble_utg_align",
        method="wfmash",
        minimap2="minimap2",
        biopython_mode="globalxx",
        mm2_threads=1,
        ps_match=2,
        ps_mismatch=-3,
        ps_open=5,
        ps_extend=2,
    ):
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
            return {"fasta": fasta, "result": paf, "cmd": cmd}

        elif method == "wfmash":
            paf = f"{output_prefix}.paf"
            cmd = f"wfmash -t {mm2_threads} -n 1 -X -N {fasta} {fasta} > {paf}"
            import subprocess

            subprocess.run(cmd, shell=True, check=True)
            return {"fasta": fasta, "result": paf, "cmd": cmd}
        elif method == "edlib":
            try:
                import edlib
            except ImportError as e:
                raise ImportError("Please intall edlib: pip install edlib") from e

            data = self.extract_bubble_internal_unique_sequences(bubble)
            a, b2 = internals
            s1 = data[a]["seq"]
            s2 = data[b2]["seq"]
            qlen, tlen = len(s1), len(s2)

            def parse_ext_cigar(cigar: str):
                m = x = ins = dele = 0
                num = []
                for ch in cigar or "":
                    if ch.isdigit():
                        num.append(ch)
                        continue
                    n = int("".join(num)) if num else 1
                    num = []
                    if ch == "=":
                        m += n
                    elif ch == "X":
                        x += n
                    elif ch == "I":
                        ins += n
                    elif ch == "D":
                        dele += n
                    elif ch == "M":
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
                    n = int("".join(num)) if num else 1
                    num = []
                    op = "M" if ch in ("=", "X", "M") else ch
                    if op == last_op:
                        run += n
                    else:
                        if last_op is not None:
                            flush(last_op, run)
                        last_op, run = op, n
                if last_op is not None:
                    flush(last_op, run)
                return "".join(out)

            def run_one(target: str, strand: str):
                res = edlib.align(s1, target, mode="HW", task="path")
                cigar = res.get("cigar") or ""
                # edlib.locations: [(start, end)] 0-based, inclusive
                locs = res.get("locations") or []
                if locs and isinstance(locs[0], (list, tuple)) and len(locs[0]) == 2:
                    ts_incl, te_incl = int(locs[0][0]), int(locs[0][1])
                else:
                    ts_incl, te_incl = 0, len(target) - 1
                m, x, ins, dele, alen, q_consume, t_consume = parse_ext_cigar(cigar)
                nmatch = m
                ident = (nmatch / alen) if alen > 0 else 0.0

                qstart, qend = 0, q_consume
                if strand == "+":
                    tstart, tend = ts_incl, te_incl + 1
                else:
                    tstart = tlen - 1 - te_incl
                    tend = (tlen - 1 - ts_incl) + 1
                return {
                    "identity": ident,
                    "nmatch": nmatch,
                    "alen": alen,
                    "qstart": qstart,
                    "qend": qend,
                    "tstart": tstart,
                    "tend": tend,
                    "strand": strand,
                    "cigar_ext": cigar,
                }

            info_f = run_one(s2, "+")
            info_r = run_one(self._revcomp2(s2), "-")
            best = (
                info_f
                if (info_f["identity"], info_f["alen"])
                >= (info_r["identity"], info_r["alen"])
                else info_r
            )

            paf = f"{output_prefix}.paf"
            mapq = 255
            cg = ext_to_cg(best["cigar_ext"]) if best["cigar_ext"] else ""
            with open(paf, "w") as fh:
                fields = [
                    a,
                    str(qlen),
                    str(best["qstart"]),
                    str(best["qend"]),
                    best["strand"],
                    b2,
                    str(tlen),
                    str(best["tstart"]),
                    str(best["tend"]),
                    str(best["nmatch"]),
                    str(best["alen"]),
                    str(mapq),
                ]
                line = "\t".join(fields)
                if cg:
                    line += f"\tcg:Z:{cg}\tNM:i:{best['alen'] - best['nmatch']}"
                fh.write(line + "\n")

            return {"fasta": fasta, "result": paf, "cmd": "edlib:HW"}

        elif method == "biopython":
            try:
                from Bio import pairwise2
            except ImportError:
                raise ImportError("Can not import biopython for method=biopython")
            seqs = {}
            for line in open(fasta):
                if line.startswith(">"):
                    cur = line[1:].split()[0]
                    seqs[cur] = []
                else:
                    seqs[cur].append(line.strip())
            seqs = {k: "".join(v) for k, v in seqs.items()}
            a, b = internals
            s1 = seqs[a]
            s2 = seqs[b]
            if biopython_mode == "globalxx":
                alns = pairwise2.align.globalxx(s1, s2, one_alignment_only=True)
            else:
                alns = pairwise2.align.localxx(s1, s2, one_alignment_only=True)
            out_txt = f"{output_prefix}.biopy.aln.txt"
            with open(out_txt, "w") as fh:
                if alns:
                    aln = alns[0]
                    fh.write(f"{a} vs {b}\nScore={aln.score}\n{aln.seqA}\n{aln.seqB}\n")
                else:
                    fh.write("No alignment produced\n")
            return {"fasta": fasta, "result": out_txt}
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
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 12:
                    continue
                qn, qlen, qs, qe, strand, tn, tlen, ts, te, nmatch, alen, mapq = parts[
                    :12
                ]
                if qn == tn:
                    continue
                if not ({qn, tn} == {a, b}):
                    continue
                qlen = _to_int(qlen)
                qs, qe = _to_int(qs), _to_int(qe)
                tlen = _to_int(tlen)
                ts, te = _to_int(ts), _to_int(te)
                nmatch = _to_int(nmatch)
                alen = _to_int(alen, default=max(qe - qs, te - ts))
                mapq = _to_int(mapq, default=0)
                ident = (nmatch / alen) if alen > 0 else 0.0
                cg = None
                for tag in parts[12:]:
                    if tag.startswith("cg:Z:"):
                        cg = tag[5:]
                        break
                rec = {
                    "qname": qn,
                    "tname": tn,
                    "strand": strand,
                    "qlen": qlen,
                    "qstart": qs,
                    "qend": qe,
                    "tlen": tlen,
                    "tstart": ts,
                    "tend": te,
                    "nmatch": nmatch,
                    "alen": alen,
                    "mapq": mapq,
                    "identity": ident,
                    "cigar": cg,
                }
                if (
                    (best is None)
                    or (rec["identity"] > best["identity"])
                    or (
                        rec["identity"] == best["identity"]
                        and rec["alen"] > best["alen"]
                    )
                ):
                    best = rec
        return best

    def _parse_edlib_cigar(self, cigar: str):
        m = x = ins = dele = 0
        num = []
        for ch in cigar:
            if ch.isdigit():
                num.append(ch)
                continue
            n = int("".join(num)) if num else 1
            num = []
            if ch == "=":
                m += n
            elif ch == "X":
                x += n
            elif ch == "I":
                ins += n
            elif ch == "D":
                dele += n
            elif ch == "M":
                m += n
        alen = m + x + ins + dele
        return m, x, ins, dele, alen

    def align_bubbles_and_collect(
        self,
        bubbles,
        out_dir="bubble_align_all",
        utg_prefix="utg",
        exact_path_len=3,
        require_degree=False,
        method="wfmash",
        minimap2="minimap2",
        biopython_mode="globalxx",
        tsv_path=None,
        n_jobs=20,
        mm2_threads=1,
    ):
        targets = self.filter_bubbles_two_internal_utg_pattern(
            bubbles,
            utg_prefix=utg_prefix,
            exact_path_len=exact_path_len,
            require_degree=require_degree,
        )
        self._ensure_dir(out_dir)

        def do_one(idx, b):
            bid = f"b{idx}"
            internals = b["internal_nodes"]
            data = self.extract_bubble_internal_unique_sequences(b)
            a, b2 = internals[0], internals[1]
            len_a = len(data[a]["seq"])
            len_b = len(data[b2]["seq"])
            meta_a = data[a]["meta"]
            meta_b = data[b2]["meta"]
            prefix = os.path.join(out_dir, f"{bid}_{b['source']}_{b['sink']}")
            try:
                res = self.align_internal_utgs(
                    b,
                    output_prefix=prefix,
                    method=method,
                    minimap2=minimap2,
                    biopython_mode=biopython_mode,
                    mm2_threads=mm2_threads,
                )
                fasta_path = res["fasta"]
                out_path = res["result"]
                if method == "minimap2":
                    paf_best = self._parse_pair_paf(out_path, a, b2)
                    if paf_best is None:
                        row = {
                            "bubble_id": bid,
                            "source": b["source"],
                            "sink": b["sink"],
                            "int1": a,
                            "int2": b2,
                            "len1": len_a,
                            "len2": len_b,
                            "leftOv1": meta_a["left_overlap"],
                            "rightOv1": meta_a["right_overlap"],
                            "leftOv2": meta_b["left_overlap"],
                            "rightOv2": meta_b["right_overlap"],
                            "method": method,
                            "fasta": fasta_path,
                            "paf": out_path,
                            "identity": np.nan,
                            "alen": 0,
                            "nmatch": 0,
                            "mapq": np.nan,
                            "strand": None,
                            "qstart": None,
                            "qend": None,
                            "tstart": None,
                            "tend": None,
                            "cigar": None,
                        }
                    else:
                        row = {
                            "bubble_id": bid,
                            "source": b["source"],
                            "sink": b["sink"],
                            "int1": a,
                            "int2": b2,
                            "len1": len_a,
                            "len2": len_b,
                            "leftOv1": meta_a["left_overlap"],
                            "rightOv1": meta_a["right_overlap"],
                            "leftOv2": meta_b["left_overlap"],
                            "rightOv2": meta_b["right_overlap"],
                            "method": method,
                            "fasta": fasta_path,
                            "paf": out_path,
                            "identity": paf_best["identity"],
                            "alen": paf_best["alen"],
                            "nmatch": paf_best["nmatch"],
                            "mapq": paf_best["mapq"],
                            "strand": paf_best["strand"],
                            "qstart": paf_best["qstart"],
                            "qend": paf_best["qend"],
                            "tstart": paf_best["tstart"],
                            "tend": paf_best["tend"],
                            "cigar": paf_best["cigar"],
                        }
                elif method in ("edlib", "parasail"):
                    st = res.get("stats", {})
                    row = {
                        "bubble_id": bid,
                        "source": b["source"],
                        "sink": b["sink"],
                        "int1": a,
                        "int2": b2,
                        "len1": len_a,
                        "len2": len_b,
                        "leftOv1": meta_a["left_overlap"],
                        "rightOv1": meta_a["right_overlap"],
                        "leftOv2": meta_b["left_overlap"],
                        "rightOv2": meta_b["right_overlap"],
                        "method": method,
                        "fasta": fasta_path,
                        "paf_or_txt": out_path,
                        "identity": st.get("identity", np.nan),
                        "alen": st.get("alen", 0),
                        "nmatch": st.get("nmatch", 0),
                        "mapq": np.nan,
                        "strand": st.get("strand", "+"),
                        "qstart": st.get("qstart"),
                        "qend": st.get("qend"),
                        "tstart": st.get("tstart"),
                        "tend": st.get("tend"),
                        "cigar": st.get("cigar"),
                    }
                else:
                    row = {
                        "bubble_id": bid,
                        "source": b["source"],
                        "sink": b["sink"],
                        "int1": a,
                        "int2": b2,
                        "len1": len_a,
                        "len2": len_b,
                        "leftOv1": meta_a["left_overlap"],
                        "rightOv1": meta_a["right_overlap"],
                        "leftOv2": meta_b["left_overlap"],
                        "rightOv2": meta_b["right_overlap"],
                        "method": method,
                        "fasta": fasta_path,
                        "aln_txt": out_path,
                    }
            except Exception as e:
                logger.warning(f"bubble {bid} failed align: {e}")
                row = {
                    "bubble_id": bid,
                    "source": b["source"],
                    "sink": b["sink"],
                    "int1": a,
                    "int2": b2,
                    "len1": len_a,
                    "len2": len_b,
                    "leftOv1": meta_a["left_overlap"],
                    "rightOv1": meta_a["right_overlap"],
                    "leftOv2": meta_b["left_overlap"],
                    "rightOv2": meta_b["right_overlap"],
                    "method": method,
                    "fasta": None,
                    "paf": None,
                    "identity": np.nan,
                    "alen": 0,
                    "nmatch": 0,
                    "mapq": np.nan,
                    "strand": None,
                    "qstart": None,
                    "qend": None,
                    "tstart": None,
                    "tend": None,
                    "cigar": None,
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
            df.to_csv(tsv_path, sep="\t", index=False)
        return df

    def write(self, out_handle):
        for name, ln in self.segments.items():
            has_ln = False
            tags = self.segment_tags.get(name, [])
            new_tags = []
            for t in tags:
                if t.startswith("LN:i:"):
                    new_tags.append(f"LN:i:{ln}")
                    has_ln = True
                else:
                    new_tags.append(t)
            if not has_ln:
                new_tags.insert(0, f"LN:i:{ln}")
            
            tags_str = "\t".join(new_tags)
            out_handle.write(f"S\t{name}\t*\t{tags_str}\n")
          
        for e in self.links:
            length1 = self.segments[e["from"]]
            length2 = self.segments[e["to"]]
            out_handle.write(
                f"L\t{e['from']}\t{e['from_orient']}\t{e['to']}\t{e['to_orient']}\t{e['overlap']}\tL1:i:{length1}\tL2:i:{length2}\n"
            )


    def generate_mock_hic_files(
        self,
        clm_path="mock.clm",
        contacts_path="mock.split.contacts",
        weight_scale=0.1,  
        min_weight=50,   
    ):
        with xopen(clm_path, "w") as f_clm, xopen(contacts_path, "w") as f_cnt:
            for e in self.links:
                ctg1 = e["from"]
                o1 = e["from_orient"]
                ctg2 = e["to"]
                o2 = e["to_orient"]

                ctg1_ht = f"{ctg1}_1" if o1 == "+" else f"{ctg1}_0"
                ctg2_ht = f"{ctg2}_0" if o2 == "+" else f"{ctg2}_1"

                weight = max(min_weight, int(e["ov_len"] * weight_scale))

                f_cnt.write(f"{ctg1_ht}\t{ctg2_ht}\t{weight}\n")

                for dir1 in ["+", "-"]:
                    for dir2 in ["+", "-"]:
                        if dir1 == o1 and dir2 == o2:
                            dist = 10
                        else:
                            dist = 100000 

                        distances = " ".join([str(dist)] * weight)

                        ctg1_dir = f"{ctg1}{dir1}"
                        ctg2_dir = f"{ctg2}{dir2}"

                        if ctg1_dir > ctg2_dir:
                            ctg1_dir, ctg2_dir = ctg2_dir, ctg1_dir

                        f_clm.write(f"{ctg1_dir} {ctg2_dir}\t{weight}\t{distances}\n")

        logger.info(
            f"Generated mock Hi-C files: `{clm_path}`, `{contacts_path}` from {len(self.links)} overlaps."
        )
        return clm_path, contacts_path


    def apply_chimeric_bed(self, bed_file):
        """
        Update the GFA graph using a BED file containing chimeric breakpoints.
        BED format: original_contig  start  end  new_contig_name
        """
        breaks = defaultdict(list)
        with open(bed_file, 'r') as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.strip().split()
                orig = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                new_name = parts[3]
                breaks[orig].append({"start": start, "end": end, "new_name": new_name})

        for orig in breaks:
            breaks[orig].sort(key=lambda x: x["start"])

        new_segments = {}
        new_seqs = {}
        new_tags = {}
        for orig, pieces in breaks.items():
            if orig in self.segments:
                orig_tags = self.segment_tags.get(orig, [])
                for p in pieces:
                    new_name = p["new_name"]
                    length = p["end"] - p["start"] + 1
                    new_segments[new_name] = length
                    new_tags[new_name] = orig_tags[:] 

                    
                    if orig in self.seqment_seqs:
                        new_seqs[new_name] = self.seqment_seqs[orig][p["start"] - 1 : p["end"]]


                self.segments.pop(orig, None)
                self.seqment_seqs.pop(orig, None)

        self.segments.update(new_segments)
        self.seqment_seqs.update(new_seqs)
        if hasattr(self, "segment_tags"):
            self.segment_tags.update(new_tags)
        new_links = []
        for link in self.links:
            frm = link["from"]
            to = link["to"]
            fo = link["from_orient"]
            to_o = link["to_orient"]

            new_frm = frm
            new_to = to

            if frm in breaks:
                if fo == "+":
                    new_frm = breaks[frm][-1]["new_name"]  # + strand out, attaches to 3' end
                else:
                    new_frm = breaks[frm][0]["new_name"]   # - strand out, attaches to 5' end
            
            if to in breaks:
                if to_o == "+":
                    new_to = breaks[to][0]["new_name"]     # + strand in, attaches to 5' end
                else:
                    new_to = breaks[to][-1]["new_name"]    # - strand in, attaches to 3' end

            link["from"] = new_frm
            link["to"] = new_to
            new_links.append(link)

        self.links = new_links
        self._graph = None  # Invalidate cached graph
        logger.info(f"Applied chimeric breaks for {len(breaks)} contigs.")


    def _get_directed_outgoing_links(self, u, u_strand):
        res = set()
        for e in self.links:
            if e["from"] == u and e["from_orient"] == u_strand:
                res.add((e["to"], e["to_orient"], e["ov_len"]))
            elif e["to"] == u:
                rev_to = "-" if e["to_orient"] == "+" else "+"
                if rev_to == u_strand:
                    rev_from = "-" if e["from_orient"] == "+" else "+"
                    res.add((e["from"], rev_from, e["ov_len"]))
        return res

    def rescue_agp_with_gfa(
        self,
        in_agp,
        out_agp,
        default_gap=100,
        rescue=True,
        remove_redundant=True,
        redundant_perc=100.0,
    ):
        from collections import defaultdict

        redundant_set = set()
        if remove_redundant:
            logger.info(
                f"Scanning for redundant contigs (Overlap >= {redundant_perc}%)..."
            )
            res_df = self.calculate_overlap_percentage()
            if not res_df.empty:
                redundant_set = set(
                    res_df[res_df["Overlap_Percentage"] >= redundant_perc]["Chromosome"]
                )
            logger.info(
                f"Identified {len(redundant_set)} redundant contigs to be removed."
            )

        agp_tigs = set()
        placed_tigs = set() 
        
        raw_scaffolds = defaultdict(list)
        with open(in_agp, "r") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.strip().split("\t")
                chr_name = p[0]
                part_type = p[4]

                if part_type in ("W", "D", "F"):
                    tig = p[5]
                    agp_tigs.add(tig)
                    is_redundant = tig in redundant_set
                    if (tig != chr_name) and (not is_redundant):
                        placed_tigs.add(tig)
                    raw_scaffolds[chr_name].append(
                        {
                            "type": "W",
                            "id": tig,
                            "tig_s": int(p[6]),
                            "tig_e": int(p[7]),
                            "strand": p[8],
                            "ov_after": 0,
                            "deleted": (tig in redundant_set),
                        }
                    )
                elif part_type in ("U", "N"):
                    raw_scaffolds[chr_name].append(
                        {"type": "GAP", "len": int(p[5]), "deleted": False}
                    )
        
        unplaced = agp_tigs - placed_tigs
        if remove_redundant and redundant_set:
            unplaced -= redundant_set

        logger.info(
            f"AGP contigs: total={len(agp_tigs)}, placed={len(placed_tigs)}, candidate_unplaced={len(unplaced)}"
        )


        out_edges = defaultdict(set)
        for e in self.links:
            frm, d_from = e["from"], e["from_orient"]
            to, d_to = e["to"], e["to_orient"]
            ov = e["ov_len"]

            out_edges[(frm, d_from)].add((to, d_to, ov))

            rev_from = "-" if d_from == "+" else "+"
            rev_to = "-" if d_to == "+" else "+"
            out_edges[(to, rev_to)].add((frm, rev_from, ov))

        scaffolds = defaultdict(list)
        corrections, overlap_trims = 0, 0
        rescued_internal = 0
        rescued_count = 0

        for chr_name, cseq in raw_scaffolds.items():
            w_nodes = [x for x in cseq if x["type"] == "W" and not x.get("deleted")]
            has_gap = any(x["type"] == "GAP" for x in cseq)
            if len(w_nodes) >= 1 or has_gap:
                scaffolds[chr_name] = cseq
                for x in w_nodes:
                    unplaced.discard(x["id"])
            else:
                continue

        logger.info(
            f"Loaded {len(scaffolds)} main scaffolds and identified {len(unplaced)} available unplaced contigs."
        )
    
        if rescue:

            for chr_name, cseq in scaffolds.items():
                w_idx = [
                    i
                    for i, x in enumerate(cseq)
                    if x["type"] == "W" and not x.get("deleted")
                ]
                insertions = {}
                for k in range(len(w_idx) - 1):
                    i1, i2 = w_idx[k], w_idx[k + 1]
                    u, u_str = cseq[i1]["id"], cseq[i1]["strand"]
                    v, v_str = cseq[i2]["id"], cseq[i2]["strand"]

                    links = list(out_edges.get((u, u_str), set()))

                    candidate_matches = [l for l in links if l[0] == v]
                    if candidate_matches:
                        candidate_matches.sort(key=lambda x: x[2], reverse=True)
                        match = candidate_matches[0]
                    else:
                        match = None


                
                    if match:
                        v_expected_str, ov = match[1], match[2]

                        if v_str != v_expected_str:
                            logger.debug(
                                f"[Correction] Adjusted orientation {chr_name}:{v} to {v_expected_str}"
                            )
                            cseq[i2]["strand"] = v_expected_str
                            corrections += 1

                        cseq[i1]["ov_after"] = ov
                        if ov > 0:
                            overlap_trims += 1
                        for g_i in range(i1 + 1, i2):
                            cseq[g_i]["deleted"] = True
                    else:
                        best_bridge = None
                        best_score = -1

                        for l_u in links:
                            u_node, u_in_strand, ov_a_u = l_u
                            if u_node in unplaced:
                                u_links = out_edges.get((u_node, u_in_strand), set())
                                # match_u = next((l for l in u_links if l[0] == v), None)

                                # if match_u:
                                #     v_expected_str, ov_u_b = match_u[1], match_u[2]

                                #     logger.info(
                                #         f"[Rescue] Filled internal GAP in {chr_name} between {u} and {v} with {u_node}"
                                #     )
                                #     rescued_internal += 1
                                #     unplaced.remove(u_node)

                                #     if v_str != v_expected_str:
                                #         logger.info(
                                #             f"[Correction] Adjusted orientation {chr_name}:{v} to {v_expected_str}"
                                #         )
                                #         cseq[i2]["strand"] = v_expected_str
                                #         corrections += 1

                                #     cseq[i1]["ov_after"] = ov_a_u
                                #     if ov_a_u > 0:
                                #         overlap_trims += 1

                                #     new_u = {
                                #         "type": "W",
                                #         "id": u_node,
                                #         "tig_s": 1,
                                #         "tig_e": self.segments[u_node],
                                #         "strand": u_in_strand,
                                #         "ov_after": ov_u_b,
                                #         "deleted": False,
                                #     }
                                #     insertions[i1] = [new_u]
                                #     if ov_u_b > 0:
                                #         overlap_trims += 1

                                #     for g_i in range(i1 + 1, i2):
                                #         cseq[g_i]["deleted"] = True
                                #     break
                                cands_u = [l for l in u_links if l[0] == v]
                                if cands_u:
                                    cands_u.sort(key=lambda x: x[2], reverse=True)
                                    match_u = cands_u[0]
                                    v_expected_str, ov_u_b = match_u[1], match_u[2]
                                    
                                    score = ov_a_u + ov_u_b
                                    if score > best_score:
                                        best_score = score
                                        best_bridge = (u_node, u_in_strand, ov_a_u, v_expected_str, ov_u_b)


                        if best_bridge:
                            u_node, u_in_strand, ov_a_u, v_expected_str, ov_u_b = best_bridge

                            logger.debug(
                                f"[Rescue] Filled internal GAP in {chr_name} between {u} and {v} with {u_node}"
                            )
                            rescued_internal += 1
                            # unplaced.remove(u_node)
                            if u_node not in redundant_set:
                                unplaced.discard(u_node)

                            if v_str != v_expected_str:
                                logger.debug(
                                    f"[Correction] Adjusted orientation {chr_name}:{v} to {v_expected_str}"
                                )
                                cseq[i2]["strand"] = v_expected_str
                                corrections += 1

                            cseq[i1]["ov_after"] = ov_a_u
                            if ov_a_u > 0:
                                overlap_trims += 1

                            new_u = {
                                "type": "W",
                                "id": u_node,
                                "tig_s": 1,
                                "tig_e": self.segments[u_node],
                                "strand": u_in_strand,
                                "ov_after": ov_u_b,
                                "deleted": False,
                            }
                            insertions[i1] = [new_u]
                            if ov_u_b > 0:
                                overlap_trims += 1

                            for g_i in range(i1 + 1, i2):
                                cseq[g_i]["deleted"] = True     

                if insertions:
                    new_cseq = []
                    for i, item in enumerate(cseq):
                        new_cseq.append(item)
                        if i in insertions:
                            new_cseq.extend(insertions[i])
                    scaffolds[chr_name] = new_cseq
            logger.info(
                f"Finished internal bridge rescue round! Bridged {rescued_internal} internal gaps."
            )

            for chr_name, cseq in scaffolds.items():
                cseq = [x for x in cseq if not x["deleted"]]

                while False:
                    last_w = next((x for x in reversed(cseq) if x["type"] == "W"), None)
                    if not last_w:
                        break

                    links = [
                        l
                        for l in out_edges.get((last_w["id"], last_w["strand"]), set())
                        if l[0] in unplaced
                    ]
                    if links:
                        links.sort(key=lambda x: x[2], reverse=True)
                        if len(links) == 1 or links[0][2] >= 2 * links[1][2]:
                            v, v_str, ov = links[0]
                            unplaced.remove(v)
                            last_w["ov_after"] = ov
                            if ov > 0:
                                overlap_trims += 1
                            cseq.append(
                                {
                                    "type": "W",
                                    "id": v,
                                    "tig_s": 1,
                                    "tig_e": self.segments[v],
                                    "strand": v_str,
                                    "ov_after": 0,
                                    "deleted": False,
                                }
                            )
                            logger.info(
                                f"[Rescue] Appended {v}({v_str}) to Tail of {chr_name}"
                            )
                            rescued_count += 1
                        else:
                            break
                    else:
                        break

                while False:
                    first_w = next((x for x in cseq if x["type"] == "W"), None)
                    if not first_w:
                        break

                    rev_s = "-" if first_w["strand"] == "+" else "+"
                    links = [
                        l
                        for l in out_edges.get((first_w["id"], rev_s), set())
                        if l[0] in unplaced
                    ]
                    if links:
                        links.sort(key=lambda x: x[2], reverse=True)
                        if len(links) == 1 or links[0][2] >= 2 * links[1][2]:
                            u, u_rev_s, ov = links[0]
                            unplaced.remove(u)
                            u_str = "-" if u_rev_s == "+" else "+"

                            new_node = {
                                "type": "W",
                                "id": u,
                                "tig_s": 1,
                                "tig_e": self.segments[u],
                                "strand": u_str,
                                "ov_after": ov,
                                "deleted": False,
                            }
                            if ov > 0:
                                overlap_trims += 1
                            cseq.insert(0, new_node)
                            logger.info(
                                f"[Rescue] Prepended {u}({u_str}) to Head of {chr_name}"
                            )
                            rescued_count += 1
                        else:
                            break
                    else:
                        break
                scaffolds[chr_name] = cseq

        if remove_redundant:
            old_unplaced_len = len(unplaced)
            unplaced = unplaced - redundant_set
            logger.info(
                f"Cleaned {old_unplaced_len - len(unplaced)} unplaced redundant contigs."
            )

            for chr_name, cseq in scaffolds.items():
                if (
                    len(cseq) == 1
                    and cseq[0]["type"] == "W"
                    and cseq[0]["id"] in redundant_set
                ):
                    cseq[0]["deleted"] = True

        u_idx = 1
        for u_tig in sorted(list(unplaced)):
            scaffolds[f"{u_tig}"] = [
                {
                    "type": "W",
                    "id": u_tig,
                    "tig_s": 1,
                    "tig_e": self.segments[u_tig],
                    "strand": "+",
                    "ov_after": 0,
                    "deleted": False,
                }
            ]
            u_idx += 1

        logger.info(
            f"Done graph traversing. Outputting AGP: Rescued={rescued_count}, Corrected={corrections}, OverlapsTrimmed={overlap_trims}."
        )

        is_handle = hasattr(out_agp, "write")
        out = out_agp if is_handle else open(out_agp, "w")
        try:
            for chr_name, cseq in scaffolds.items():
                cleaned_cseq = []
                for x in cseq:
                    if x.get("deleted"):
                        continue
                    if x["type"] == "GAP":
                        if not cleaned_cseq or cleaned_cseq[-1]["type"] == "GAP":
                            continue
                    cleaned_cseq.append(x)

                if cleaned_cseq and cleaned_cseq[-1]["type"] == "GAP":
                    cleaned_cseq.pop()

                if not cleaned_cseq:
                    continue
                obj_s = 1
                part = 1
                for i, item in enumerate(cleaned_cseq):
                    if item.get("deleted"):
                        continue
                    if item["type"] == "GAP":
                        gap_len = item["len"]
                        obj_e = obj_s + gap_len - 1
                        out.write(
                            f"{chr_name}\t{obj_s}\t{obj_e}\t{part}\tU\t{gap_len}\tscaffold\tyes\tmap\n"
                        )
                        obj_s = obj_e + 1
                        part += 1
                    elif item["type"] == "W":
                        tig_s, tig_e = item["tig_s"], item["tig_e"]

                        if i > 0:
                            prev = cleaned_cseq[i - 1]
                            if prev["type"] == "W" and prev.get("ov_after", 0) > 0:
                                O = prev["ov_after"]
                                if item["strand"] == "+":
                                    tig_s += O
                                else:
                                    tig_e -= O

                        if tig_s > tig_e:
                            tig_s, tig_e = 1, 1

                        c_len = tig_e - tig_s + 1
                        obj_e = obj_s + c_len - 1

                        out.write(
                            f"{chr_name}\t{obj_s}\t{obj_e}\t{part}\tW\t{item['id']}\t{tig_s}\t{tig_e}\t{item['strand']}\n"
                        )
                        obj_s = obj_e + 1
                        part += 1
        finally:
            if not is_handle:
                out.close()

        out_name = getattr(out_agp, "name", str(out_agp)) if is_handle else str(out_agp)
        logger.info(f"Rescued AGP exported to `{out_name}`")

        if out_name != "<stdout>":
            if is_handle:
                out_agp.flush()
            from .cli import statagp
            try:
                statagp.main(args=[out_name, "-o", f"{out_name}.stat"], prog_name="statagp")
            except SystemExit as e:
                exc_info = sys.exc_info()
                exit_code = e.code
                if exit_code is None:
                    exit_code = 0
                
                if exit_code != 0:
                    raise e

        return out_agp

    
    def order_and_orient_group(self, contig_group):
        group = set(contig_group)
        edges = defaultdict(list)
        in_degree = {n: 0 for n in group}

        for e in self.links:
            u, u_str = e["from"], e["from_orient"]
            v, v_str = e["to"], e["to_orient"]
            ov = e["ov_len"]

            if u in group and v in group:
                edges[(u, u_str)].append((v, v_str, ov))
                in_degree[v] += 1
                
                rev_u = "-" if u_str == "+" else "+"
                rev_v = "-" if v_str == "+" else "+"
                edges[(v, rev_v)].append((u, rev_u, ov))

        for k in edges:
            edges[k].sort(key=lambda x: x[2], reverse=True)

        all_paths = []
        visited = set()


        starts = [n for n in group if in_degree[n] == 0]
        if not starts:
            starts = sorted(list(group), key=lambda x: self.segments.get(x, 0), reverse=True)
        else:
            starts = sorted(starts, key=lambda x: self.segments.get(x, 0), reverse=True)

        for start_node in starts:
            if start_node in visited:
                continue

            curr = start_node
            curr_str = "+"
            if len(edges.get((curr, "-"), [])) > len(edges.get((curr, "+"), [])):
                curr_str = "-"

            current_path = []
            while curr and curr not in visited:
                current_path.append((curr, curr_str))
                visited.add(curr)
            
                next_node, next_str = None, None
                for nxt, nxt_str, ov in edges.get((curr, curr_str), []):
                    if nxt not in visited:
                        next_node = nxt
                        next_str = nxt_str
                        break
                
                curr, curr_str = next_node, next_str
            
            if current_path:
                all_paths.append(current_path)

        for node in group:
            if node not in visited:
                all_paths.append([(node, "+")])
                visited.add(node)

        return all_paths