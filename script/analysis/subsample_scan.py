import pandas as pd
import numpy as np
import os
import itertools
import argparse
from scipy.special import comb
from time import strftime
import glob

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--conditiona", help="RNA-seq data for Condition A",
                    type=str, required=True)
parser.add_argument("-b", "--conditionb", help="RNA-seq data for Condition B",
                    type=str, required=True)
parser.add_argument("-o", "--output", help="output directory", type=str,
                    default='{}_subsample_output'.format(strftime('%Y%m%d')))
parser.add_argument("-prefix", "--prefix",
                    help="Add prefix to the log",
                    type=str, default='')

def get_combinations(df1, df2, n1, n2, counter=0, times=50):
    s1 = df1.shape[1]
    s2 = df2.shape[1]
    combinations = comb(s1, n1)*comb(s2, n2)
    # check whether to do all the combinations or do random subsampling
    if combinations < 50:
        comb1 = list(itertools.combinations(df1.columns, n1))
        comb2 = list(itertools.combinations(df2.columns, n2))
        pairs = list(itertools.product(comb1, comb2))
        # path = '/home/shared/Data/epee/data/rnaseq/coad_subsampling'
        for i, p in enumerate(pairs):
            out_df1 = df1.loc[:, p[0]]
            out_df2 = df2.loc[:, p[1]]
            path1 = '{}/{}_subsample_{}_cond1.txt'.format(outdir, counter,
                                                          args.prefix)
            path2 = '{}/{}_subsample_{}_cond2.txt'.format(outdir, counter,
                                                          args.prefix)
            path1_exists = os.path.exists(path1)
            path2_exists = os.path.exists(path2)
            if path1_exists and path2_exists:
                continue
            else:
                out_df1.to_csv(path1, sep='\t')
                out_df2.to_csv(path2, sep='\t')
            counter += 1
    else:
        for i in range(times):
            np.random.seed(i)
            path1 = '{}/{}_subsample_{}_cond1.txt'.format(outdir, counter,
                                                          args.prefix)
            path2 = '{}/{}_subsample_{}_cond2.txt'.format(outdir, counter,
                                                          args.prefix)
            path1_exists = os.path.exists(path1)
            path2_exists = os.path.exists(path2)
            if path1_exists and path2_exists:
                continue
            else:
                out_df1 = df1.loc[:, np.random.choice(df1.columns, n1)]
                out_df2 = df2.loc[:, np.random.choice(df2.columns, n2)]
                out_df1.to_csv(path1, sep='\t')
                out_df2.to_csv(path2, sep='\t')
            counter += 1
    return counter


if __name__ == '__main__':
    args = parser.parse_args()

    outdir = '{d}{p}_subsampling/'.format(d=args.output, p=args.prefix)
    # Read files to subsample
    df1 = pd.read_csv(args.conditiona, sep='\t', index_col=0)
    df2 = pd.read_csv(args.conditionb, sep='\t', index_col=0)

    # Generate subsample files if not generated already
    os.makedirs(os.path.dirname(outdir), exist_ok=True)

    # Th2
    if args.prefix == 'Th2':
        c = 0
        for n in range(1, 6):
            c = get_combinations(df1, df2, n, n, counter=c)

    # COAD
    if args.prefix == 'COAD':
        _ = get_combinations(df1, df2,
                             int(df1.shape[1]*0.7),
                             int(df2.shape[1]*0.7))

    # Write sbatch script to run the commands in GPU partition
