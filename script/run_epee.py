# Copyright 2018, Murat Can Cobanoglu
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# filter future warnings
import warnings
warnings.simplefilter("ignore", category=FutureWarning)

from epee import *
import numpy as np
import pandas as pd
import argparse
import logging
import time
import os
import itertools
import multiprocessing
from time import localtime, strftime

# set tensorflow verbosity
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'



parser = argparse.ArgumentParser()
parser.add_argument("-a", "--conditiona", help="RNA-seq data for Condition A",
                    type=str, required=True)
parser.add_argument("-b", "--conditionb", help="RNA-seq data for Condition B",
                    type=str, required=True)
parser.add_argument("-na", "--networka", help="Network for condition A",
                    type=str, required=True)
parser.add_argument("-nb", "--networkb", help="Network for condition B",
                    type=str, required=True)
# DEFAULTS
parser.add_argument("-o", "--output", help="output directory", type=str,
                    default='')
parser.add_argument("-reg1", "--lregularization", help="lasso regularization \
                    parameter", type=float, default=0.01)
parser.add_argument("-reg2", "--gregularization", help="graph contrained \
                    regularization parameter", type=float, default=0.01)
parser.add_argument("-s", "--step", help="optimizer learning-rate",
                    type=float, default=0.0001)
parser.add_argument("-c", "--conditioning", help="Weight for the interactions \
                    not known", type=bool, default=True)
parser.add_argument("-r", "--runs", help="Number of independent runs", type=int,
                    default=20)
parser.add_argument("-i", "--iterations", help="Number of iterations",
                    type=int, default=100000)
parser.add_argument("-norm", "--normalize", help="""
                    Weight normalization strategy. Default:"minmax"
                    Valid options: {"minmax", "log", "log10", "no"} """,
                    type=str, default='minmax')
parser.add_argument("-model", "--model", help="""
                    Model regularization choice. Default: "epee-gcl"
                    Valid options: {"epee-gcl","epee-l","no-penalty" """,
                    type=str, default='epee-gcl')
parser.add_argument("-verbose", "--verbose",
                    help="logging info levels 10, 20, or 30",
                    type=int, default=10)
# OPTIONAL SETTINGS
parser.add_argument("-eval", "--evaluate",
                    help="Evaluation mode available for Th1, Th2, Th17, \
                    Bmem, COAD, and AML",
                    type=str, default=None)
parser.add_argument("-prefix", "--prefix",
                    help="Add prefix to the log",
                    type=str, default=strftime('%Y%m%d'))
# OPTIONAL FLAGS
parser.add_argument("-weight", "--weight",
                    help="Store all the inferred weights",
                    action='store_true')
parser.add_argument("-multiprocess", "--multiprocess",
                    help="multiprocess the calculation of perturb and \
                    regulator scores", action='store_true')
# NULL FLAG
parser.add_argument("-null", "--null",
                    help="Generate null scores by label permutation",
                    action='store_true')
# NULL SETTINGS
parser.add_argument("-seed", "--seed", help="Starting seed number",
                    type=int, default=0)
parser.add_argument("-perturb", "--perturb", help="True label perturb scores. Required when running permutations for null model",
                    type=str, default=None)


def get_scores(sel):
    """To get perturb and regulator score"""
    y1, w1, w1_df, y2, w2, w2_df, count = sel

    # Calculate perturb scores
    genescore_runi = get_perturb_scores(Y1, y1, X1, w1,
                                        Y2, y2, X2, w2, S1, S2)
    genescore_runi.columns = ['gene', 'set{}'.format(count)]

    if args.null:
        regscore_runi, diff_regs = get_diff_regulatory_activity(
                                         actual_perturb['gene'][:1000],
                                         w1_df, w2_df, top_regs=20)
    else:
        regscore_runi, diff_regs = get_diff_regulatory_activity(
                                         genescore_runi['gene'][:1000],
                                         w1_df, w2_df, top_regs=20)

    regscore_runi.columns = ['gene', 'set{}'.format(count)]

    return (genescore_runi, regscore_runi)


def run_epee():
    """To run EPEE with specified inputs."""
    logging.info('SAMPLES: Y1: {} | Y2: {}'.format(Y1.shape[1], Y2.shape[1]))
    logging.info('Tensorflow: {}'.format(tf.__version__))
    logging.info('GENES: {}'.format(Y1.shape[0]))
    logging.info('TFs: {}'.format(S1.shape[1]))
    logging.info('MODEL LEARNING STARTED')
    genescore_df = pd.DataFrame()
    regscore_df = pd.DataFrame()
    loss_runs = []
    y1_s = []
    y2_s = []
    w1_s = []
    w2_s = []
    w1S1_s = []
    w2S2_s = []

    for rid in range(args.runs):
        start = time.time()
        logging.debug('Tensorflow: {}'.format(tf.__version__))
        logging.debug('MODEL: {} learning Y1'.format(rid))
        y1, w1, loss_arr1 = run_model(np.array(Y1), np.array(X1),
                                      np.array(S1),
                                      l_reg=args.lregularization,
                                      g_reg=args.gregularization,
                                      step=args.step,
                                      itr=args.iterations,
                                      log_itr=round(args.iterations/20),
                                      seed=rid+args.seed,
                                      model=args.model,
                                      val=condition_val)
        logging.debug('MODEL: {} learning Y2'.format(rid))
        y2, w2, loss_arr2 = run_model(np.array(Y2), np.array(X2),
                                      np.array(S2),
                                      l_reg=args.lregularization,
                                      g_reg=args.gregularization,
                                      step=args.step,
                                      itr=args.iterations,
                                      log_itr=round(args.iterations/20),
                                      seed=rid+args.seed,
                                      model=args.model,
                                      val=condition_val)

        loss_runs.append((rid, loss_arr1[-1], loss_arr2[-1]))

        # Calculate w1S1 and w2S2
        w1_s1 = np.multiply(w1, S1)
        w2_s2 = np.multiply(w2, S2)

        w1_df = get_weights_df(w1_s1, Y1.index, X1.index)
        w2_df = get_weights_df(w2_s2, Y2.index, X2.index)

        w1o_df = get_weights_df(w1, Y1.index, X1.index)
        w2o_df = get_weights_df(w2, Y2.index, X2.index)

        # Store dataframes
        y1_s.append(y1)
        y2_s.append(y2)
        w1_s.append(w1)
        w2_s.append(w2)
        w1S1_s.append(w1_df)
        w2S2_s.append(w2_df)

        # Output inferred weights if args.weight is True and args.null is False
        if args.weight and not args.null:
            w1o_df.to_csv('{}/model/w1_{}.txt'.format(outdir, rid),
                          sep='\t')
            w2o_df.to_csv('{}/model/w2_{}.txt'.format(outdir, rid),
                          sep='\t')
            if rid == 0:
                S1.to_csv('{}/model/S1_input.txt'.format(outdir),
                          sep='\t')
                S2.to_csv('{}/model/S2_input.txt'.format(outdir),
                          sep='\t')
                X1.to_csv('{}/model/X1_input.txt'.format(outdir),
                          sep='\t')
                X2.to_csv('{}/model/X2_input.txt'.format(outdir),
                          sep='\t')
                Y1.to_csv('{}/model/Y1_input.txt'.format(outdir),
                          sep='\t')
                Y2.to_csv('{}/model/Y2_input.txt'.format(outdir),
                          sep='\t')

        end = time.time()

        logging.info('MODEL: {} RUNTIME: {} mins'.format(rid,
                     round((end-start)/60, 2)))

    # For each pairs of inferred weights calculate perturb and regulator scores
    # logging.info('CALCULATE PERTURB AND REGULATOR SCORES')
    logging.info('SCORES: pairwise comparision of all Y1 and Y2 models')

    list_runs = list(range(args.runs))
    pairs = list(itertools.product(list_runs, list_runs))
    score_inputs = []
    for count, p in enumerate(pairs):
        m1, m2 = p
        score_inputs.append((y1_s[m1], w1_s[m1], w1S1_s[m1],
                             y2_s[m2], w2_s[m2], w2S2_s[m2],
                             count))

    if args.multiprocess:
        cpu_count = multiprocessing.cpu_count()
        p = multiprocessing.Pool(int(cpu_count/2))
        out = p.map(get_scores, score_inputs)
    else:
        out = []
        for i in score_inputs:
            i_out = get_scores(i)
            out.append(i_out)

    for count, scores in enumerate(out):
        genescore_runi, regscore_runi = scores
        if count == 0:
            genescore_df = genescore_runi.copy()
            regscore_df = regscore_runi.copy()
        else:
            # if np.all(genescore_runi.index == genescore_df.index):
            #     genescore_df[genescore_runi.columns[1]] = genescore_runi.iloc[:, 1]
            # else:
            genescore_df = pd.merge(genescore_df, genescore_runi, on='gene')
            # if np.all(regscore_runi.index == regscore_df.index):
            #     regscore_df[regscore_runi.columns[1]] = regscore_runi.iloc[:, 1]
            # else:
            regscore_df = pd.merge(regscore_df, regscore_runi, on='gene')

    sum_genescore_df = get_summary_scoresdf(genescore_df)
    sum_regscore_df = get_summary_scoresdf(regscore_df)

    if args.null:
        sum_regscore_df.to_csv('{}/null/regulator_scores_{}.txt'.format(
                               outdir, args.seed),
                               sep='\t')
        sum_genescore_df.to_csv('{}/null/perturb_scores_{}.txt'.format(
                                outdir, args.seed),
                                sep='\t')
        regscore_df.to_csv('{}/null/all_regulator_scores_{}.txt'.format(
                           outdir, args.seed),
                           sep='\t')
        genescore_df.to_csv('{}/null/all_perturb_scores_{}.txt'.format(
                            outdir, args.seed),
                            sep='\t')
    else:
        sum_regscore_df.to_csv('{}/scores/regulator_scores.txt'.format(
                               outdir), sep='\t')
        sum_genescore_df.to_csv('{}/scores/perturb_scores.txt'.format(
                                outdir), sep='\t')
        regscore_df.to_csv('{}/scores/all_regulator_scores.txt'.format(
                           outdir), sep='\t')
        genescore_df.to_csv('{}/scores/all_perturb_scores.txt'.format(
                            outdir), sep='\t')

        loss_df = pd.DataFrame(loss_runs)
        loss1_df = pd.DataFrame(loss_arr1)
        loss2_df = pd.DataFrame(loss_arr2)
        loss_df.to_csv('{}/model/loss_runs.txt'.format(outdir),
                       sep='\t')
        loss1_df.to_csv('{}/model/loss1_arr_y1.txt'.format(outdir),
                        sep='\t')
        loss2_df.to_csv('{}/model/loss2_arr_y2.txt'.format(outdir),
                        sep='\t')

    if args.evaluate:

        known_regs = {'Th2': ['STAT6', 'GATA3'],
                      'Th1': ['TBX21', 'STAT1'],
                      'Th17': ['ARID5A', 'RORA', 'STAT3'],
                      'Bmem': ['STAT5A'],
                      'AML': ['WT1', 'MYB', 'ETV6', 'SOX4', 'CEBPA', 'RUNX1'],
                      'COAD': ['MYC', 'KLF4']
                      }
        regs = X1.index
        actual = known_regs[args.evaluate]
        line = [args.model, args.prefix, args.seed, args.lregularization,
                args.gregularization, len(Y1.columns), len(Y2.columns),
                args.conditiona, args.conditionb]

        for t in actual:
            line.append(list(regscore_df['gene']).index(t))

        logging.info(line)

        median_rank = get_median_rank_position(actual,
                                               list(regscore_df['gene']),
                                               len(regs), args.model)
        logging.info('RANK: {}'.format(median_rank[0]))
        print('RANK: {}'.format(median_rank[0]))

        line.append(median_rank[0])

        performance = '\t'.join(map(str, line))


        outfilename = '{}/{}_eval.txt'.format(args.output, args.evaluate)
        # write header if filename does not exists
        if not os.path.isfile(outfilename):
            with open(outfilename, 'a') as myfile:
                headerlist = ['MODEL', 'PREFIX', 'SEED', 'L1', 'L2',
                              'SamplesY1', 'SamplesY2', 'PathY1', 'PathY2']
                for t in actual:
                    headerlist.append('{}_RANK'.format(t))
                headerlist.append('MEDIAN_RANK')
                header = '\t'.join(headerlist)
                myfile.write('{}\n'.format(header))

        with open(outfilename, 'a') as myfile:
            myfile.write('{}\n'.format(performance))


if __name__ == '__main__':

    run_start = time.time()
    args = parser.parse_args()
    outdir = '{d}{p}_{m}_{l}_{g}_{s}'.format(d=args.output, p=args.prefix,
                                                  m=args.model,
                                                  l=args.lregularization,
                                                  g=args.gregularization,
                                                  s=args.seed)
    if args.null:
        if args.perturb == None:
            raise RuntimeError('Please provide perturb genes output generated with actual labels. --perturb <path to perturb_scores.txt>')
        os.makedirs(os.path.dirname('{}/null/'.format(outdir)),
                    exist_ok=True)
        logging.basicConfig(filename='{}/null_log.txt'.format(outdir),
                            level=args.verbose)
        actual_perturb = pd.read_csv(args.perturb, index_col=0, sep='\t')
    else:
        os.makedirs(os.path.dirname('{}/model/'.format(outdir)),
                    exist_ok=True)
        os.makedirs(os.path.dirname('{}/scores/'.format(outdir)),
                    exist_ok=True)
        logging.basicConfig(filename='{}/log.txt'.format(outdir),
                            level=args.verbose)

    logging.info('####### {} STARTING ANALYSIS  #######'.format(args.prefix))
    logging.debug('Multiprocessing: {}'.format(args.multiprocess))
    logging.info('EPEE: {}'.format(strftime("%a, %d %b %Y %H:%M:%S",
                 localtime())))
    Y1, Y2, X1, X2, S1, S2, condition_val = get_epee_inputs(
                                     args.conditiona, args.conditionb,
                                     args.networka, args.networkb,
                                     conditioning=args.conditioning,
                                     weightNormalize=args.normalize,
                                     null=args.null, seed=args.seed)
    run_epee()
    run_end = time.time()

    logging.info('Time elapsed: {} mins'.format(
                                            round((run_end-run_start)/60, 2)))
    logging.info('#######  {} ANALYSIS COMPLETED ########'.format(args.prefix))
