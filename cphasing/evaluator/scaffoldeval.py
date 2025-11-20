#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
evaluation of the contig orientation and ordering based on simulation data
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np

from ..agp import import_agp
from collections import defaultdict



def get_adjacencies(paths):
    adjacencies = set()
    for _, path in paths.items():
        if len(path) < 2:
            continue
        for i in range(len(path) - 1):
            c1 = path[i]
            c2 = path[i + 1]
            id1, strand1 = c1[:-1], c1[-1]
            id2, strand2 = c2[:-1], c2[-1]
            node1 = f"{id1}_tail" if strand1 == "+" else f"{id1}_head"
            node2 = f"{id2}_head" if strand2 == "+" else f"{id2}_tail"
            adj = tuple(sorted((node1, node2)))
            adjacencies.add(adj)
    return adjacencies

def get_id_adjacencies(paths):
    s = set()
    for _, path in paths.items():
        if len(path) < 2:
            continue
        for i in range(len(path) - 1):
            c1 = path[i]
            c2 = path[i + 1]
            id1 = c1[:-1]
            id2 = c2[:-1]
            pair = tuple(sorted((id1, id2)))
            s.add(pair)
    return s

def get_oriented_adjacencies_by_id(paths):
  
    m = {}
    for _, path in paths.items():
        if len(path) < 2:
            continue
        for i in range(len(path) - 1):
            c1 = path[i]
            c2 = path[i + 1]
            id1, s1 = c1[:-1], c1[-1]
            id2, s2 = c2[:-1], c2[-1]
            n1 = f"{id1}_tail" if s1 == "+" else f"{id1}_head"
            n2 = f"{id2}_head" if s2 == "+" else f"{id2}_tail"
            key = tuple(sorted((id1, id2)))
            oriented = tuple(sorted((n1, n2)))
            m[key] = oriented
    return m

def evaluate_adjacencies(test_adj, truth_adj):
    tp_set = test_adj.intersection(truth_adj)
    fp_set = test_adj.difference(truth_adj)
    fn_set = truth_adj.difference(test_adj)
    tp = len(tp_set)
    fp = len(fp_set)
    fn = len(fn_set)
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1_score = (
        2 * (precision * recall) / (precision + recall)
        if (precision + recall) > 0
        else 0.0
    )
    return {
        "true_positives": tp,
        "false_positives (misjoins)": fp,
        "false_negatives (breaks)": fn,
        "precision": precision,
        "recall": recall,
        "f1_score": f1_score,
        "misjoins_list": sorted(list(fp_set)),
    }





def _normalize_group_name(s: str) -> str:
    if s is None:
        return ""
    x = str(s).strip()
    x = x.replace("|", "_")
    x_low = x.lower()
    for pref in (
        "chr",
        "chromosome_",
        "chromosome",
        "scaffold_",
        "scaf_",
        "ctg_",
        "contig_",
    ):
        if x_low.startswith(pref):
            x = x[len(pref) :]
            x_low = x_low[len(pref) :]
            break
    # 去掉前导 0
    x = x.lstrip("0")
    # 常见同义
    if x_low in ("mt", "m", "mitochondria"):
        return "MT"
    if x_low in ("x", "y"):
        return x.upper()
    return x


def _build_contig_group_map(paths, normalize=True):
    """
    paths: dict[group_name -> [contig_with_strand,...]]
    return: dict[contig_id -> group_name]
    """
    m = {}
    for chrom, lst in paths.items():
        g = _normalize_group_name(chrom) if normalize else chrom
        for c in lst:
            cid = c[:-1]  
            m[cid] = g
    return m


def infer_group_alias_by_overlap(truth_paths, test_paths, normalize=True):
    tmap = _build_contig_group_map(truth_paths, normalize=normalize)
    pmap = _build_contig_group_map(test_paths, normalize=normalize)

    counts = {}  # truth_g -> {test_g: n}
    common = set(tmap.keys()) & set(pmap.keys())
    for cid in common:
        tg = tmap[cid]
        pg = pmap[cid]
        d = counts.setdefault(tg, {})
        d[pg] = d.get(pg, 0) + 1

    alias = {}
    for tg, d in counts.items():
        pg = max(d.items(), key=lambda kv: (kv[1], kv[0]))[0]
        alias[tg] = pg
    return alias


def compute_grouping_errors(truth_paths, test_paths, alias=None, normalize=True):
    tmap = _build_contig_group_map(truth_paths, normalize=normalize)
    pmap = _build_contig_group_map(test_paths, normalize=normalize)

    if alias is None:
        alias = infer_group_alias_by_overlap(
            truth_paths, test_paths, normalize=normalize
        )

    mapped_tmap = {cid: alias.get(g, g) for cid, g in tmap.items()}

    common = set(mapped_tmap.keys()) & set(pmap.keys())
    if not common:
        return 0, 0.0
    errors = sum(1 for cid in common if mapped_tmap[cid] != pmap[cid])
    rate = errors / len(common) if common else 0.0
    return errors, rate


def _build_contig_lengths(*dfs):
    lens = {}
    for df in dfs:
        if df is None or "type" not in df.columns:
            continue
        for _, row in df[df["type"] == "W"].iterrows():
            cid = row["id"]
            L = None
            if "size" in row and pd.notnull(row["size"]):
                L = int(row["size"])
            elif "length" in row and pd.notnull(row["length"]):
                L = int(row["length"])
            elif (
                {"start", "end"}.issubset(row.index)
                and pd.notnull(row["start"])
                and pd.notnull(row["end"])
            ):
                L = int(row["end"]) - int(row["start"]) + 1
            if L is not None:
                if cid not in lens or L > lens[cid]:
                    lens[cid] = L
    return lens


def get_weighted_adjacencies(paths, contig_len):
    adj_w = {}
    for _, path in paths.items():
        if len(path) < 2:
            continue
        for i in range(len(path) - 1):
            c1 = path[i]
            c2 = path[i + 1]
            id1, s1 = c1[:-1], c1[-1]
            id2, s2 = c2[:-1], c2[-1]
            n1 = f"{id1}_tail" if s1 == "+" else f"{id1}_head"
            n2 = f"{id2}_head" if s2 == "+" else f"{id2}_tail"
            adj = tuple(sorted((n1, n2)))
            w = contig_len.get(id1, 0) + contig_len.get(id2, 0)
            adj_w[adj] = max(adj_w.get(adj, 0), w)
    return adj_w


def evaluate_adjacencies_weighted(test_adj_w, truth_adj_w):
    test_keys = set(test_adj_w.keys())
    truth_keys = set(truth_adj_w.keys())
    inter = test_keys & truth_keys
    sum_pred = float(sum(test_adj_w.values()))
    sum_truth = float(sum(truth_adj_w.values()))
    tp_w_pred = float(sum(test_adj_w[k] for k in inter)) if inter else 0.0
    tp_w_truth = float(sum(truth_adj_w[k] for k in inter)) if inter else 0.0
    w_precision = tp_w_pred / sum_pred if sum_pred > 0 else 0.0
    w_recall = tp_w_truth / sum_truth if sum_truth > 0 else 0.0
    w_f1 = (
        (2 * w_precision * w_recall / (w_precision + w_recall))
        if (w_precision + w_recall) > 0
        else 0.0
    )
    return {"w_precision": w_precision, "w_recall": w_recall, "w_f1": w_f1}


def scaf_eval(truth_agp, test_agp):
    truth_agp_df, _ = import_agp(truth_agp, split=True)
    test_agp_df, _ = import_agp(test_agp, split=True)
    
    truth_agp_df.reset_index(inplace=True)
    test_agp_df.reset_index(inplace=True)
    test_agp_df = test_agp_df[test_agp_df['chrom'] != test_agp_df['id']]

    truth_paths = defaultdict(list)
    test_paths = defaultdict(list)

    grouped_contigs = set()
    for _, row in test_agp_df.iterrows():
        if row["type"] == "W":
            contig_str = f"{row['id']}{'+' if row['orientation'] == '+' else '-'}"
            test_paths[row["chrom"]].append(contig_str)
            grouped_contigs.add(row["id"])

    for _, row in truth_agp_df.iterrows():
        if row["type"] == "W":
            if row["id"] not in grouped_contigs:
                continue
            contig_str = f"{row['id']}{'+' if row['orientation'] == '+' else '-'}"
            truth_paths[row["chrom"]].append(contig_str)
    

    truth_adjacencies = get_adjacencies(truth_paths)
    test_adjacencies = get_adjacencies(test_paths)
    results = evaluate_adjacencies(test_adjacencies, truth_adjacencies)
    group_err_cnt, group_err_rate = compute_grouping_errors(truth_paths, test_paths)

    lens = _build_contig_lengths(truth_agp_df, test_agp_df)
    res_w = evaluate_adjacencies_weighted(
        get_weighted_adjacencies(test_paths, lens),
        get_weighted_adjacencies(truth_paths, lens),
    )
    id_adj_truth = get_id_adjacencies(truth_paths)
    id_adj_test = get_id_adjacencies(test_paths)
    order_fp = sorted(id_adj_test - id_adj_truth)
    order_fn = sorted(id_adj_truth - id_adj_test)
    id_tp = id_adj_test & id_adj_truth

    ori_truth = get_oriented_adjacencies_by_id(truth_paths)
    ori_test = get_oriented_adjacencies_by_id(test_paths)
    orient_err = [p for p in id_tp if ori_test.get(p) != ori_truth.get(p)]

    print(
        "TP\tFP\tFN\tPrecision\tRecall\tF1\tMisjoinsCount\tGroupErrors\tGroupErrRate\twPrecision\twRecall\twF1"
    )

    print(
        "\t".join(
            [
                str(results["true_positives"]),
                str(results["false_positives (misjoins)"]),
                str(results["false_negatives (breaks)"]),
                f"{results['precision']:.4f}",
                f"{results['recall']:.4f}",
                f"{results['f1_score']:.4f}",
                str(len(results.get("misjoins_list", []))),
                str(group_err_cnt),
                f"{group_err_rate:.4f}",
                f"{res_w.get('w_precision', 0.0):.4f}",
                f"{res_w.get('w_recall', 0.0):.4f}",
                f"{res_w.get('w_f1', 0.0):.4f}",
            ]
        )
    )

    print("OrderFP\tOrderFN\tOrientErr")
    print(f"{len(order_fp)}\t{len(order_fn)}\t{len(orient_err)}")


    # output order and orientation errors
    # if results.get("misjoins_list", []):
    #     print("\nMisjoins (false positives):")
    #     for adj in results["misjoins_list"]:
    #         print(f"{adj[0]}\t{adj[1]}")
    # if order_fp:
    #     print("\nOrder false positives:")
    #     for pair in order_fp:
    #         print(f"{pair[0]}\t{pair[1]}")
    if order_fn:
        print("\nOrder errors:")
        for pair in order_fn:
            print(f"{pair[0]}\t{pair[1]}", file=sys.stderr)
    


    if orient_err:
        print("\nOrientation errors:")
        for pair in orient_err:
            # o1 = ori_test.get(pair)
            # o2 = ori_truth.get(pair)
            print(f"{pair[0]}\t{pair[1]}", file=sys.stderr)
    
