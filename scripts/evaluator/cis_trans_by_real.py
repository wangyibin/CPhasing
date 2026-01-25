#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
calculate the cis/trans ratio according from real order cool
"""

import argparse
import logging
import os
import os.path as op
import sys

import cooler 
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np 

from itertools import combinations, groupby
from math import sqrt
from statannotations.Annotator import Annotator

def wi_test(data1, data2):
    """
    Wilcoxon rank-sum tests
    return: pvalue
    """
    from scipy import stats
    wi = stats.ranksums(data1, data2)
    return wi.pvalue

def as_si(x, ndp):
    """
    Scientific count format function, use x replace e.
    """
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pOpt.add_argument('-c', '--cool', 
            nargs="+",
            help='Path to cool file',
            required=True)
    pOpt.add_argument("-hl", "--homo_list",
                    
            help="homolog chroms rellationship",
            required=True)
    pOpt.add_argument("--labels", nargs="*", default=None,
                      help="labels of different cool data, the number of labels must equal to cool")
    pOpt.add_argument("--color_index", nargs="*", default=None, 
                      help="color index used for each labels")
    pOpt.add_argument("--group", nargs="*", default=None,
                      help="group different input sample")
    pOpt.add_argument('--wi_test', action="store_true", default=False,
                      help="plot the wi-test, only support two cool file")
    pOpt.add_argument("--rotation", default=0, type=int, 
                      help="rotate the xticks [default: %(default)s]")
    pOpt.add_argument('-o', '--output', type=str,
            default="boxplot.png", help='output picture [default: "boxplot.png"]')
    pOpt.add_argument('-h', '--help', action='help',
        help='show help message and exit.')
    
    args = p.parse_args(args)

    if args.labels:
        assert len(args.cool) == len(args.labels), "number of labels must equal to cool files"

    results = []

    for i, cool_file in enumerate(args.cool):
        cool = cooler.Cooler(cool_file)
        matrix = cool.matrix(balance=False, sparse=True)

        homo_chroms = [i.strip().split() for i in open(args.homo_list) if i.strip()]

        total_cis = 0 
        total_trans = 0
        res = []
        pair_num = 0
        for homo in homo_chroms:
            homo = sorted(homo)
            for pair in combinations(homo, 2):
                pair_num += 1
                try:
                    cis1 = matrix.fetch(pair[0]).sum()
                    cis2 = matrix.fetch(pair[1]).sum()
                    trans = matrix.fetch(*pair).sum()
                except:
                    continue 
                    
                # data = trans / sqrt(cis1 * cis2)
                data = trans / (trans + cis1 + cis2)


                res.append(data)

        results.append(res)

    res_df = pd.DataFrame(results).T
    if args.labels:
        res_df.columns = args.labels
    

    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype'] = 42
    fig, ax = plt.subplots(figsize=(5, 7))
    boxprops_r = dict(facecolor='#a83836',color='black', linewidth=1.5)
    boxprops_b = dict(facecolor='#265e8a',color='black', linewidth=1.5)
    boxprops = dict(color='black', linewidth=1.5)
    medianprops=dict(color='black', linewidth=2.5)
    whiskerprops = dict(linestyle='--')

    
    # a = ax.boxplot(results[0], showfliers=False, patch_artist=True, notch=True, widths=0.4,
    #            boxprops=boxprops_r, medianprops=medianprops, whiskerprops=whiskerprops)

    # a_upper_extreme = [ item.get_ydata() for item in a['whiskers']][1][1]
    # a_bottom_extreme = [ item.get_ydata() for item in a['whiskers']][0][1]

    # b = ax.boxplot([[], results[1]], showfliers=False, patch_artist=True, notch=True, widths=0.4,
    #            boxprops=boxprops_b, medianprops=medianprops, whiskerprops=whiskerprops)
    # b_upper_extreme = [ item.get_ydata() for item in b['whiskers']][3][1]
    # b_bottom_extreme = [ item.get_ydata() for item in b['whiskers']][2][1]

    # # for patch, color in zip(bplot['boxes'], ['#a83836',  '#8dc0ed', '#df8384', ]):
    # #     patch.set_facecolor(color)
    # max_upper = max(a_upper_extreme, b_upper_extreme)
    # min_upper = min(a_upper_extreme, b_upper_extreme)
    # min_bottom = min(a_bottom_extreme, b_bottom_extreme)
    # h = max_upper/5

    # if max_upper == a_upper_extreme:
    #     y1 = max_upper + max_upper/10
    #     y2 = min_upper + max_upper/10
    # else:
    #     y1 = min_upper + max_upper/10
    #     y2 = max_upper + max_upper/10

    # ax.plot([1,1,2,2], [y1, max_upper + h, max_upper + h, y2], linewidth=1.0, color='black')
    # ax.text((1 + 2)*.5, max_upper + max_upper/4.5, r'$P = {0:s}$'.format(as_si(pvalue, 2)), ha='center', va='bottom' )
    colors = ['#df8384', '#8dc0ed', '#a83836',  "#253761",]
    # colors = ['#a83836', '#a83836', '#a83836', '#a83836',  
    #             '#df8384', '#df8384', '#df8384', '#df8384', 
    #             "#253761", '#8dc0ed',]
    # colors = ['#df8384',  '#8dc0ed',]
    # colors = ['#df8384']

    pairs = list(combinations(args.labels, 2))
    
    res_df = res_df.melt()
    if args.group:
        res_df['group'] = np.nan
        for group, label in zip(args.group, args.labels):
            
            res_df.loc[res_df['variable'] == label, 'group'] = group 
   
    if args.color_index:
        new_color = []
        for i in args.color_index:
            new_color.append(colors[int(i)])
    else:
        new_color = colors

    if args.group:
        hue_order = 'group'
    else:
        hue_order = None

    res_df.to_csv(args.output.replace(".png", "") + ".csv", sep='\t', index=True, header=True)
    print(res_df)
    ax = sns.violinplot(data=res_df, x='variable',
                        y='value', hue_order=hue_order, ax=ax, 
                        palette=new_color, alpha=1, legend=True)
    if args.wi_test:
        annotator = Annotator(ax, pairs, data=res_df, x='variable', y='value')
        annotator.configure(test="Mann-Whitney", text_format='full', show_test_name=False, loc='inside')
        annotator.apply_and_annotate()

    ax.set_xticks(range(len(args.cool)))
    # ax.set_xticklabels(['MAPQ>=1', "MAPQ>=2"])
    # ax.set_xticklabels(['Raw', 'Removed \n$\mathit{trans}$'])
    ax.set_xticklabels(args.labels, fontsize=8, rotation=args.rotation, ha='right')
    ax.set_xlabel("")
    max_y = max(max(results)) * 1.5
    # plt.ylim(-0.001, max_y)
    # ax.set_xticklabels(["Before\n realign", "After\n realign", ], fontsize=20)
    
    # ax.set_xticks([0, 1, 2])
    # ax.set_xticklabels(["Before\n realign", "After\n realign", "After\n realign2"], fontsize=20)

    # ax.set_xticks([0, 1, 2, 3, 4])
    # ax.set_xticklabels([2, 4, 6, 8, 12], fontsize=20)

    # ax.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    # ax.set_xticklabels(["k15 w5", "k15_w10", "k17_w7", "k27_w14", "k15 w5", "k15_w10", "k17_w7", "k27_w14", "Hi-C $\it{Dpn}$II", "Hi-C Arima"], 
                    #    fontsize=8, rotation=45, ha="right")
    # ax.set_xticks([0, 1, 2, 3])
    # ax.set_xticklabels(["Pore-C\n$\it{Hin}$dIII", "Pore-C\n$\it{Dpn}$II", "Hi-C\n$\it{Dpn}$II", "Hi-C\nArima"], fontsize=20)

    # ax.set_xticks([0, 1])
    # ax.set_xticklabels(["Pore-C", "Hi-C", ], fontsize=20)
    # max_y = max(max(results)) * 1.3
    # plt.ylim(-0.1, max_y)
 
    

    # ax.set_xticks([0, 1, 2, 3])
    # ax.set_xticklabels(["Pore-C\n$\it{Hin}$dIII", "Pore-C\n$\it{Dpn}$II", "Hi-C\n$\it{Dpn}$II", "Hi-C\nArima"], fontsize=20)
    ylim = ax.get_ylim()
    if args.group:
        import matplotlib.patches as mpatches
        
        unique_groups = []
        seen = set()
        group_color_map = {}
        
        for grp, col in zip(args.group, new_color):
            if grp not in seen:
                unique_groups.append(grp)
                seen.add(grp)
                group_color_map[grp] = col
        
        legend_handles = []
        for grp in unique_groups:
            patch = mpatches.Patch(facecolor=group_color_map[grp], 
                                   label=grp, edgecolor='black',
                                   linewidth=1.5)
            legend_handles.append(patch)

        ax.legend(handles=legend_handles, title="", loc='best', frameon=False, fontsize=12)
    sns.despine()
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    plt.tick_params(which='both', width=1.5, length=5)

    plt.text(0, ylim[1] * 0.05, f"n={pair_num}", fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=16)
    ax.set_ylabel(r"h-$\mathit{trans}$ error rate", fontsize=20)
    plt.savefig(args.output, dpi=600, bbox_inches='tight')
    plt.savefig(args.output.replace("png", "pdf"), dpi=600, bbox_inches='tight')

if __name__ == "__main__":
    main(sys.argv[1:])
