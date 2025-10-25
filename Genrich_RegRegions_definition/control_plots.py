#/usr/bin/env/python

"""
Script producing customized quality control plots for the crossregmap workflow.
"""

import argparse
import os

import numpy as np

import seaborn as sns
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

sns.set_palette("muted")

def load_bed_region_size(input_file):
    """
    Loads a bed file in a python dict.
    key = col 4 (region name) and a values = dict with chr, start, stop, length entries
    """
    d = {}
    with open(input_file, "r") as infile:
        for i, line in enumerate(infile):
            line_split = line.strip().split('\t')
            if len(line_split) == 1:
                line_split = line.strip().split()
            line = line_split
            ch, start, stop = line[0:3]
            d[i] = {'chr': ch, 'start': int(start), 'stop': int(stop), 'lg': int(stop) - int(start)}
    return d


def barplot_regnb(dataset, title='', out=None, gb=""):

    data_dict, lg_dict = {}, {}
    cl_dict = {}
    if gb:
        assert len(gb.keys()) == 1
        key = list(gb.keys())[0]
    for data, name in dataset:
        if gb:
            for val in gb[key]:
                if val in name:
                    cl_dict[val] = cl_dict.get(val, [])
                    if name not in cl_dict[val]:
                        cl_dict[val].append(name)
                    break
        data_dict[name] = load_bed_region_size(data)
        lg_dict[name] = len(data_dict[name])

    df = pd.DataFrame.from_dict(lg_dict, orient='index', columns=['Number of peaks'])
    df['Sample'] = df.index

    if gb:
        ORDER = []
        COLORS = []
        PALETTE = sns.color_palette("muted").as_hex()
        for i, val in enumerate(sorted(list(cl_dict.keys()))):
            ORDER += cl_dict[val]
            COLORS += [PALETTE[i] for _ in cl_dict[val]]

        if "mapped" not in list(gb.values())[0]:
            ax = sns.barplot(x="Sample", data=df, y='Number of peaks', order=ORDER, palette=COLORS)
        else:
            ax = sns.barplot(x="Sample", data=df, y='Number of peaks', palette=COLORS)
    else:
        ax = sns.barplot(x="Sample", data=df, y='Number of peaks')

    for p in ax.patches:
        try:
            ax.annotate(int(p.get_height()), (p.get_x() + p.get_width() / 2., p.get_height()),
                        ha='center', va='center', fontsize=8, color='gray', xytext=(0, 5),
                        textcoords='offset points')
        except:
            pass

    plt.xticks(rotation=45, horizontalalignment="right")
    sns.despine()
    plt.title(title)
    plt.tight_layout()

    if out:
        plt.savefig(out)
    else:
        plt.show()

    plt.close('all')
    return data_dict, df


def boxplot_regsize(dataset, fliers=True, xlab='', ylab='', out=None, title='', gb=None):

    def df_reg_size(data, name):
        d = load_bed_region_size(data)
        rec = [(i, float(d[i]['lg']/1000)) for i in d]
        return pd.DataFrame.from_records(rec, columns=["peaks", name])

    cl_dict = {}
    for _, name in dataset:
        if gb:
            for val in gb[key]:
                if val in name:
                    cl_dict[val] = cl_dict.get(val, [])
                    if name not in cl_dict[val]:
                        cl_dict[val].append(name)
                    break

    df = df_reg_size(dataset[0][0], dataset[0][1])

    for (data, name) in dataset[1:]:
        df_tmp = df_reg_size(data, name)
        df = df.merge(df_tmp, on='peaks', how='left')

    df.set_index("peaks", inplace=True)

    if gb:
        ORDER = []
        COLORS = []
        PALETTE = sns.color_palette("muted").as_hex()
        for i, val in enumerate(sorted(list(cl_dict.keys()))):
            ORDER += cl_dict[val]
            COLORS += [PALETTE[i] for _ in cl_dict[val]]

        if "mapped" not in list(gb.values())[0]:
            sns.boxplot(data=df, showfliers=fliers, order=ORDER, palette=COLORS)
        else:
            sns.boxplot(data=df, showfliers=fliers, palette=COLORS)


    else:
        sns.boxplot(data=df, showfliers=fliers)
    plt.xticks(rotation=45, horizontalalignment="right")
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    sns.despine()
    plt.tight_layout()

    if out:
        plt.savefig(out)
    else:
        plt.show()

    plt.close('all')

def load_tss_dist(input_file):
    d = {}
    with open(input_file, "r") as infile:
        for i, line in enumerate(infile):
            line = line.strip().split('\t')
            dist = int(line[-1])
            if dist != -1:
                if dist != 0:
                    dist = np.log(dist)
                d[i] = dist
    return d

def distplot_tss_dist(dataset, xlab='', ylab='', out=None, title=''):

    def df_tss_dist(data, name):
        tss_dist_file = os.path.dirname(data) + "/TSS_"+os.path.splitext(os.path.basename(data))[0]
        if not os.path.dirname(data):
            tss_dist_file = tss_dist_file.strip('/')
        d = load_tss_dist(tss_dist_file)
        rec = list(d.items())

        return pd.DataFrame.from_records(rec, columns=["peaks", name])

    df = df_tss_dist(dataset[0][0], dataset[0][1])

    for (data, name) in dataset[1:]:
        df_tmp = df_tss_dist(data, name)
        df = df.merge(df_tmp, on='peaks', how='outer')

    df.set_index("peaks", inplace=True)

    # plt.figure(figsize=(10,7))
    for col in df.columns:
        ax = sns.kdeplot(df[col], label=col)

    ticks = [0] + list(np.logspace(1, 9, num=9-1+1, base=10, dtype='int'))
    ticks_val = [0] + [np.log(i) for i in ticks[1:]]
    plt.xticks(ticks_val)
    lab = [int(i/1000) for i in ticks]
    lab[1] = ''
    lab[2] = ''
    lab[-2] = ''
    lab[-3] = ''
    ax.set(xticklabels=lab)
    plt.xlim([0, ticks_val[-1]])
    plt.ylim([0, 0.5])
    plt.xticks(rotation=45, horizontalalignment="right")
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.legend()
    plt.title(title)
    sns.despine()
    plt.tight_layout()
    if out:
        plt.savefig(out)
    else:
        plt.show()
    plt.close('all')
    return df.median()


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-i', '--input', nargs='+', required=True)

    PARSER.add_argument('-l', '--labels', nargs='+', required=False)

    PARSER.add_argument('-plot_types', "--plot_types", nargs='+', required=False,
                        default=["bar", "box", "tss"])

    PARSER.add_argument('-o', "--output", nargs='+', required=False, default="")

    PARSER.add_argument('-tss', '--tss', type=str, required=False)

    PARSER.add_argument('-t', '--title', type=str, required=False, default='')

    PARSER.add_argument('-gb', '--groupby', type=str, required=False, default='')

    PARSER.add_argument('--verbose', action='store_true')

    ARGS = vars(PARSER.parse_args())

    if not ARGS["output"]:
        ARGS["output"] = ARGS["plot_types"]

    if not ARGS["labels"]:
        ARGS["labels"] = [os.path.splitext(os.path.basename(i))[0] for i in ARGS["input"]]

    assert len(ARGS["labels"]) == len(ARGS["input"]),\
           "Error : different numbers of inputs and labels, please check your args."

    DATA = list(zip(ARGS["input"], ARGS["labels"]))
    TO_PLOT = {plot: {"out": ARGS["output"][i], "title": ARGS["title"]}\
               for (i, plot) in enumerate(ARGS["plot_types"])}

    gb = ""
    if ARGS["groupby"]:
        gb = {}
        key, values = ARGS["groupby"].split("=")
        gb[key] = values.split(',')

    if "bar" in TO_PLOT:
        barplot_regnb(DATA, out=TO_PLOT["bar"]["out"], title=TO_PLOT["bar"]["title"], gb=gb)

    if "box" in TO_PLOT:
        boxplot_regsize(DATA, fliers=False, ylab='Region size (kb)', out=TO_PLOT["box"]["out"],
                        title=TO_PLOT["box"]["title"], gb=gb)

    if "tss" in TO_PLOT:
        for infile, name in DATA:
            outname = os.path.dirname(infile)+"/TSS_"+os.path.splitext(os.path.basename(infile))[0]
            if not os.path.dirname(infile):
                outname = outname.strip('/')
            print(outname)
            tss_file = ARGS["tss"]
            #FIXME do this in snakefile, it seems sometimes this breaks?
            cmd1 = f"cat {infile} | tr ' ' '\t' | bedtools sort > {infile}_temp"
            # os.system(cmd1)
            cmd = cmd1 + " && " +f"bedtools closest -d -a {infile}_temp -b {tss_file} -t first > {outname} && \
                    rm {infile}_temp"
            os.system(cmd)

        distplot_tss_dist(DATA, xlab='Distance to nearest TSS (kb)', out=TO_PLOT["tss"]["out"],
                          title=TO_PLOT["tss"]["title"])

        for infile, name in DATA:
            outname = os.path.dirname(infile)+"/TSS_"+os.path.splitext(os.path.basename(infile))[0]
            if os.path.exists(outname):
                os.remove(outname)

    plt.close("all")
