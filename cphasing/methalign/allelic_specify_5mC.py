#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
from dataclasses import dataclass
from typing import List, Tuple, Dict, Iterable, Optional, Set

import pandas as pd


@dataclass
class Seg:
    # 一个匹配段（来自 CIGAR 的 M/= /X），用于双向映射
    q_start: int
    q_end: int
    t_start: int
    t_end: int
    rev: bool  # 目标是否为负向（t 坐标递减）


def parse_cigar_ops(cg: str) -> List[Tuple[str, int]]:
    ops: List[Tuple[str, int]] = []
    n = 0
    for ch in cg.strip():
        if ch.isdigit():
            n = n * 10 + ord(ch) - 48
        else:
            if n == 0:
                raise ValueError(f"Malformed CIGAR around op {ch} in {cg}")
            ops.append((ch, n))
            n = 0
    if n != 0:
        raise ValueError(f"Trailing number without op in {cg}")
    return ops


def build_segments(qstart: int, qend: int, tstart: int, tend: int, strand: str, cg: str) -> List[Seg]:
    # 将 CIGAR 展开为匹配段（M/= /X），并记录 q/t 两端区间；插入/缺失只推进坐标，不产生段
    ops = parse_cigar_ops(cg)
    segs: List[Seg] = []
    q = qstart
    if strand == '+':
        t = tstart
        for op, l in ops:
            if op in ('M', '=', 'X'):
                segs.append(Seg(q, q + l, t, t + l, False))
                q += l
                t += l
            elif op == 'I':
                q += l
            elif op in ('D', 'N'):
                t += l
    else:
        t = tend
        for op, l in ops:
            if op in ('M', '=', 'X'):
                segs.append(Seg(q, q + l, t - l, t, True))
                q += l
                t -= l
            elif op == 'I':
                q += l
            elif op in ('D', 'N'):
                t -= l
    # 基本健壮性检查
    if q != qend:
        # 有些 cg 可能来自压缩后的比对（剪切不体现在 cg），此处不强求相等
        pass
    return segs


def parse_paf_row(paf_path: str, contig1: str, contig2: str) -> pd.Series:
    # 读取 PAF，筛选 contig1->contig2 的最佳一条（按 block 或 identity）
    # 兼容含/不含 header 的简单制表符文件
    df = pd.read_csv(paf_path, sep='\t', header=None, comment='#', dtype={0: str, 5: str}, engine='c')
    # 标准列（minimap2/wfmash）：0..11 固定，12+ 为 tags
    if df.shape[1] < 12:
        raise ValueError("PAF file has fewer than 12 columns")
    df = df[df[0].astype(str) == contig1]
    df = df[df[5].astype(str) == contig2]
    if df.empty:
        raise ValueError(f"No PAF rows linking {contig1} to {contig2}")

    # 提取 cg:Z:tag
    def extract_cg(row) -> Optional[str]:
        for j in range(12, len(row)):
            val = row.iloc[j]
            if isinstance(val, str) and val.startswith('cg:Z:'):
                return val[5:]
        # 兼容某些 wfmash 导出把 cg 放在第15列（你的 read_paf 示例）
        if len(row) >= 16 and isinstance(row.iloc[15], str):
            return row.iloc[15]
        return None

    df = df.copy()
    df['cg'] = df.apply(extract_cg, axis=1)
    df = df[~df['cg'].isna()]
    if df.empty:
        raise ValueError("No cg:Z: CIGAR found in PAF rows for the given contigs")

    # 计算 identity 与 block 作为优选指标
    qstart, qend = 2, 3
    tstart, tend = 7, 8
    match_col, block_col = 9, 10
    df['identity'] = df[match_col].astype(float) / df[block_col].astype(float)
    # 选 identity 最大、其后按 block 次序
    best = df.sort_values(['identity', block_col], ascending=[False, False]).iloc[0]
    return best

def parse_paf_row_any(paf_path: str, contig1: str, contig2: str) -> Tuple[pd.Series, bool]:
    """
    返回 (best_row, swapped)
    swapped=False 表示找到了 contig1->contig2；
    swapped=True  表示找到了 contig2->contig1，需要在构建段时交换 q/t。
    """
    try:
        row = parse_paf_row(paf_path, contig1, contig2)
        return row, False
    except Exception:
        row = parse_paf_row(paf_path, contig2, contig1)
        return row, True

def load_bedgraph(bg_path: str, contig: str, value_cutoff: float) -> Set[int]:
    # 读取 bedGraph：chrom start end value
    cols = ['chrom', 'start', 'end', 'value']
    bg = pd.read_csv(bg_path, sep='\t', header=None, comment='#', names=cols, dtype={0: str})
    bg = bg[bg['chrom'] == contig]
    if value_cutoff is not None:
        bg = bg[bg['value'] >= value_cutoff]
    # 仅保留 C 坐标（start）
    return set(bg['start'].astype(int).tolist())


def map_positions_q2t(pos_list: Iterable[int], segs: List[Seg], rev_shift_c_on_target: bool) -> Dict[int, int]:
    # 将 query 的一组坐标映射到 target（不可映射返回中略，不收录）
    # 为线性时间，先排序，再按段推进
    positions = sorted(set(int(p) for p in pos_list))
    out: Dict[int, int] = {}
    si = 0
    for p in positions:
        # 推进到覆盖 p 的段
        while si < len(segs) and p >= segs[si].q_end:
            si += 1
        if si >= len(segs):
            break
        s = segs[si]
        if p < s.q_start or p >= s.q_end:
            continue
        # 在段内映射
        if not s.rev:
            t = s.t_start + (p - s.q_start)
        else:
            # 段内反向：t 从 t_end-1 递减
            t = s.t_end - 1 - (p - s.q_start)
            if rev_shift_c_on_target:
                # 把 G 坐标转成 CpG 的 C 坐标
                t -= 1
                if t < 0:
                    continue
        out[p] = t
    return out


def map_positions_t2q(pos_list: Iterable[int], segs: List[Seg], rev_shift_c_on_target: bool) -> Dict[int, int]:
    # 把 target 的坐标映射回 query；逻辑与上面对称
    positions = sorted(set(int(p) for p in pos_list))
    # 为了按 t 端推进，先按 t_start 排序的段副本
    segs_t = sorted(segs, key=lambda s: (s.t_start, s.t_end))
    out: Dict[int, int] = {}
    si = 0
    for p in positions:
        pt = p
        # 若目标是反向段，输入的是 C 坐标，段内比对对应到 G，需要先把 C->G（+1）
        # 这样才能落入 [t_start, t_end)
        # 与 q2t 的 "-1" 对称
        # 注意：只对 rev 段处理，正向段不需要
        hit = False
        while si < len(segs_t) and p >= segs_t[si].t_end:
            si += 1
        k = si
        while k < len(segs_t) and p >= segs_t[k].t_start:
            s = segs_t[k]
            if s.rev:
                # 输入 p 是 C 坐标，段内匹配到的是 G 坐标 → +1
                pt = p + 1 if rev_shift_c_on_target else p
            else:
                pt = p
            if pt >= s.t_start and pt < s.t_end:
                if not s.rev:
                    q = s.q_start + (pt - s.t_start)
                else:
                    # 反向段：q 增加时 t 递减
                    q = s.q_start + (s.t_end - 1 - pt)
                out[p] = q
                hit = True
                break
            k += 1
        if not hit:
            continue
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Identify allele-specific 5mC sites between two contigs using wfmash PAF (cg) and bedGraph."
    )
    ap.add_argument("--paf", required=True, help="wfmash/minimap2 PAF with cg:Z:CIGAR")
    ap.add_argument("--bg", required=True, help="bedGraph with columns: chrom start end value")
    ap.add_argument("--contig1", required=True, help="query contig in PAF")
    ap.add_argument("--contig2", required=True, help="target contig in PAF")
    ap.add_argument("--out", required=True, help="output TSV for ASM sites")
    ap.add_argument("--value-cutoff", type=float, default=50, help="min value to keep a site from bedGraph")
    args = ap.parse_args(argv)

    # 选出 contig1->contig2 的最佳一条比对
    row, swapped = parse_paf_row_any(args.paf, args.contig1, args.contig2)
   
    strand = str(row[4])
    cg = row['cg']

    if not swapped:
        qstart, qend = int(row[2]), int(row[3])
        tstart, tend = int(row[7]), int(row[8])
        # 构建匹配段：contig1 作为 query，contig2 作为 target
        segs = build_segments(qstart, qend, tstart, tend, strand, cg)
        clip_q = (qstart, qend)
        clip_t = (tstart, tend)
    else:
        # PAF 为 contig2->contig1，交换 q/t 起止，使得 segs 仍为 contig1(作为 query) -> contig2(作为 target)
        t_qstart, t_qend = int(row[2]), int(row[3])   # 原始行的 query 区间（contig2）
        t_tstart, t_tend = int(row[7]), int(row[8])   # 原始行的 target 区间（contig1）
        # 交换后：q = 原 target(contig1)，t = 原 query(contig2)
        qstart, qend = t_tstart, t_tend
        tstart, tend = t_qstart, t_qend
        segs = build_segments(qstart, qend, tstart, tend, strand, cg)
        clip_q = (qstart, qend)
        clip_t = (tstart, tend)


    # 读甲基化位点（C 坐标）
    sites1 = load_bedgraph(args.bg, args.contig1, args.value_cutoff)
    sites2 = load_bedgraph(args.bg, args.contig2, args.value_cutoff)

    # 裁剪到对齐区间
    sites1 = set(p for p in sites1 if clip_q[0] <= p < clip_q[1])
    sites2 = set(p for p in sites2 if clip_t[0] <= p < clip_t[1])

    # contig1 -> contig2（负链时把 G->C 做 -1 平移）
    q2t = map_positions_q2t(sites1, segs, rev_shift_c_on_target=True)
    mapped1 = set(q2t.values())

    # contig2 -> contig1（负链时把 C->G 做 +1 平移，见实现内注释）
    t2q = map_positions_t2q(sites2, segs, rev_shift_c_on_target=True)
    mapped2 = set(t2q.values())

    # ASM：在一方甲基化，另一方对应位点不甲基化（且可映射）
    asm_from_1 = sorted(p for p in sites1 if p in q2t and q2t[p] not in sites2)
    asm_from_2 = sorted(p for p in sites2 if p in t2q and t2q[p] not in sites1)

    # 也可输出一致位点（可选）
    consistent = []
    for p in sites1:
        if p in q2t and q2t[p] in sites2:
            consistent.append((p, q2t[p]))

    # 输出
    rows_out = []
    for p in asm_from_1:
        rows_out.append({
            'contig': args.contig1,
            'pos': p,
            'mapped_contig': args.contig2,
            'mapped_pos': q2t[p],
            'type': 'ASM_contig1_methylated_only'
        })
    for p in asm_from_2:
        rows_out.append({
            'contig': args.contig2,
            'pos': p,
            'mapped_contig': args.contig1,
            'mapped_pos': t2q[p],
            'type': 'ASM_contig2_methylated_only'
        })
    # 可选：一致位点
    for (p1, p2) in consistent:
        rows_out.append({
            'contig': args.contig1,
            'pos': p1,
            'mapped_contig': args.contig2,
            'mapped_pos': p2,
            'type': 'consistent_both_methylated'
        })

    out_df = pd.DataFrame(rows_out,
                          columns=['contig', 'pos', 'mapped_contig', 'mapped_pos', 'type'])
    out_df.to_csv(args.out, sep='\t', index=False)

    # 简要统计打印
    print(f"#pairs: 1->{args.contig2} mapped: {len(mapped1)}, 2->{args.contig1} mapped: {len(mapped2)}")
    print(f"ASM from {args.contig1}: {len(asm_from_1)}")
    print(f"ASM from {args.contig2}: {len(asm_from_2)}")
    print(f"Consistent methylated pairs: {len(consistent)}")


if __name__ == '__main__':
    main()