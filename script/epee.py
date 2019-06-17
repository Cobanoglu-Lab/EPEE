# Copyright 2018, Viren Amin, Murat Can Cobanoglu
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import pandas as pd
import tensorflow as tf
import numpy as np
import logging
import sklearn.preprocessing as pp
import glob
import statsmodels
import statsmodels.api as sm


def _normalizeweight(df, metric):
    """To normalize the TF-gene regulatory network weights.

    Parameters
    ----------
    df
        Pandas dataframe containing network. Tab-delimited file where first
        column contains TF name, second column contains regulated gene name,
        and third column containing TF-gene regulator score.
    metric
        metric used to normalize the TF-gene regulatory score. 'minmax':
        performs the minmax scaling of regulatory scores between 1 and 2.
        'log': performs natural logarithm transformation of the regulatory
        scores. 'log10': performs base 10 logarithmic transformation of the
        regulatory score. 'no': performs no normalization to the regulatory
        score.

    Returns
    -------
    df
        Pandas dataframe containing the normalized weight

    """
    if metric == 'minmax':
        weights = np.log(df['Score']).values.reshape(-1, 1)
        norm_score = pp.MinMaxScaler(feature_range=(1, 2)).fit_transform(
                                                           weights)
    if metric == 'log':
        norm_score = np.log(df['Score']+1)
    if metric == 'log10':
        norm_score = np.log10(df['Score']+1)
    if metric == 'no':
        norm_score = df['Score']
    df['Score'] = norm_score
    return df


def _network_df(filename, Y, conditioning, weightNormalize='minmax'):
    """To take network filename and generate pandas dataframe object.

    Function takes the file path of the TF-gene regulatory network data and
    returns the pandas dataframe object containing row as gene and column as
    TF. Value in the table represents the putative regulation score of the
    TF regulating the gene.

    Parameters
    ----------
    filename
        Path to the TF-gene regulatory network data.
    Y
        Pandas dataframe containing expression data.
    conditioning
        Boolean. If conditioning is True, then TF-gene with no putative
        regulation are given small uniform score. 1/10 of the minimum TF-gene
        regulatory score.
    weightNormalize
        metric used to normalize the TF-gene regulatory score. 'minmax':
        performs the minmax scaling of regulatory scores between 1 and 2.
        'log': performs natural lograithmic transformation of the regulatory
        scores. 'log10': performs base 10 logarithmic transformation of the
        regulatory score. 'no': performs no normalization to the regulatory
        score.

    Returns
    -------
    S
        Pandas dataframe containing network. Columns are TFs and rows are
        regulated genes. Values are the weights of the TF-gene regulation.
    min_val
        value used for the conditioning

    """
    S = pd.read_csv(filename, sep='\t', header=None)
    S.columns = ['TF', 'Target', 'Score']
    S = _normalizeweight(S, weightNormalize)
    if conditioning:
        min_val = np.min(S['Score'])/10
    else:
        min_val = 0
    S = pd.pivot_table(S, values='Score', index='Target', columns='TF',
                       fill_value=min_val)
    regmask = [tf in Y.index for tf in S.columns]

    S = S.ix[Y.index, regmask].fillna(min_val)
    S = S.astype(np.float32)
    return (S, min_val)


def _exp_df(filename):
    """To take expression data filename and generate pandas dataframe object.

    Function takes the file path of the expression data and returns the
    log(normalized expression data + 1).

    Parameters
    ----------
    filename
        Path to the expression data.

    Returns
    -------
    Y
        Pandas dataframe containing log(normalized expression data + 1)

    """
    Y = pd.read_csv(filename, sep='\t', index_col=0)
    Y = np.log(Y + 1)
    Y = Y.ix[Y.index.drop_duplicates(keep=False), :]
    Y = Y.astype(np.float32)
    return Y


def _eval_indices(Y1, Y2, S1, S2):
    """To eval whether the indices are same.

    Evalulate whether indices are same for both conditions expression data and
    the context-specific network dataframes.

    Parameters
    ----------
    Y1
        Pandas dataframe containing condition 1 expression data
    Y2
        Pandas dataframe containing condition 2 expression data
    S1
        Pandas dataframe containing network 1. Columns are TFs and rows are
        regulated genes. Values are the weights of the TF-gene regulation.
    S2
        Pandas dataframe containing network 2. Columsn are TFs and rows are
        regulated genes. Values are the weights of the TF-gene regulation.

    Returns
    -------
    None

    """
    if all(Y1.index == Y2.index):
        logging.debug('Y index are equal')
    else:
        raise RuntimeError('Y1 Y2 index are not equal')

    if all(S1.index == S2.index):
        logging.debug('S index are equal')
    else:
        raise RuntimeError('S1 S2 index are not equal')

    if all(S1.index == Y1.index):
        logging.debug('S and Y are equal')
    else:
        raise RuntimeError('S Y index are not equal')


def _get_min_samples(Y1, Y2, X1, X2, min_samples, seed=0):
    """Select min samples among conditions.

    Function takes the expression dataframes of two conditions and shuffle the
    labels. If the number of labels are not evenly distrubted, then the output
    contains minimum samples of the two conditions.

    Parameters
    ----------
    Y1
        Pandas dataframe containing condition 1 expression data
    Y2
        Pandas dataframe containing condition 2 expression data
    X1
        Pandas dataframe containing TF expression across condition 1 samples
    X2
        Pandas dataframe containing TF expression across condition 2 samples
    min_samples
        Minimum number of samples between condition 1 and condition 2
        expression data
    seed
        Seed for random sampling. Default = 0.

    Returns
    -------
    sY1
        Pandas dataframe containing shuffled labels expression data for
        condition 1
    sY2
        Pandas dataframe containing shuffled labels expression data for
        condition 2
    sX1
        Pandas dataframe containing shuffled labels TF expression data
        for condition 1
    sX2
        Pandas dataframe containing shuffled labels TF expression data
        for condition 2

    """
    Y1_samples = list(Y1.columns)
    Y2_samples = list(Y2.columns)
    np.random.shuffle(Y1_samples)
    np.random.shuffle(Y2_samples)
    mY1 = Y1.loc[:, Y1_samples[:min_samples]]
    mY2 = Y2.loc[:, Y2_samples[:min_samples]]
    mX1 = X1.loc[:, Y1_samples[:min_samples]]
    mX2 = X2.loc[:, Y2_samples[:min_samples]]
    return (mY1, mY2, mX1, mX2)


def _shuffled_inputs(Y1, Y2, X1, X2, seed=0, shuffleGenes=False):
    """Shuffle the labels of the inputs samples.

    Function takes the expression dataframes of two conditions and shuffle the
    labels. 

    Parameters
    ----------
    Y1
        Pandas dataframe containing condition 1 expression data
    Y2
        Pandas dataframe containing condition 2 expression data
    X1
        Pandas dataframe containing TF expression across condition 1 samples
    X2
        Pandas dataframe containing TF expression across condition 2 samples
    seed
        Seed for random sampling. Default = 0.
    shuffleGenes
        Boolean flag for shuffling genes instead of labels. Default = False

    Returns
    -------
    sY1
        Pandas dataframe containing shuffled labels expression data for
        condition 1
    sY2
        Pandas dataframe containing shuffled labels expression data for
        condition 2
    sX1
        Pandas dataframe containing shuffled labels TF expression data
        for condition 1
    sX2
        Pandas dataframe containing shuffled labels TF expression data
        for condition 2

    """
    np.random.seed(seed)

    Y1_2 = pd.merge(Y1, Y2, left_index=True, right_index=True, how='inner')
    X1_2 = pd.merge(X1, X2, left_index=True, right_index=True, how='inner')

    n_cols = Y1.shape[1]
    if not shuffleGenes:
        # We will shuffle labels
        samples = list(Y1_2.columns)
        np.random.shuffle(samples)
        sY1 = Y1_2.loc[:, samples[:n_cols]]
        sY2 = Y1_2.loc[:, samples[n_cols:]]
        sX1 = X1_2.loc[:, samples[:n_cols]]
        sX2 = X1_2.loc[:, samples[n_cols:]]
    else:
        # We will shuffle genes
        genes = list(Y1_2.index)
        np.random.shuffle(genes)
        sY1 = Y1_2.loc[genes, Y1.columns] 
        sY2 = Y1_2.loc[genes, Y2.columns]
        for df in [sY1, sY2]:
            df.index = Y1_2.index
        tfs = list(X1_2.index)
        sX1 = sY1.loc[tfs, X1.columns]
        sX2 = sY2.loc[tfs, X2.columns]
        for df in [sX1, sX2]:
            df.index = tfs 

    return (sY1, sY2, sX1, sX2)


def get_epee_inputs(c1, c2, n1, n2, conditioning=True, weightNormalize='minmax',
                    null=False, shuffleGenes=False, seed=0):
    """To generate inputs for EPEE.

    Function takes the network and expression data filenames and generates
    inputs for EPEE to run.

    Parameters
    ----------
    c1
        filename containing condition 1 expression data
    c2
        filename containing condition 2 expression data
    n1
        filename containing network 1
    n2
        filename containing network 2
    conditioning
        whether to provide small weight value to non TF-gene putative
        interactions
    weightNormalize
        three metric of weight normalize are implemented. 'minmax' whether
        to log normalize the weights and scale the weights between 1
        and 2. 'log' whether to log normalize weight+1. 'log10' whether to
        use log base 10 normalize weight+1.
    null
        If flag is set, then samples in condition1 and condition2 are shuffled
    shuffleGenes
        If flag is set, we run null tests by shuffling the genes instead of labels
    seed
        Seed for random sampling

    Returns
    -------
    Y1
        Pandas dataframe containing condition 1 expression data
    Y2
        Pandas dataframe containing condition 2 expression data
    X1
        Pandas dataframe containing TF expression across condition 1 samples
    X2
        Pandas dataframe containing TF expression across condition 2 samples
    S1
        Pandas dataframe containing network 1. Columns are TFs and rows are
        regulated genes. Values are the weights of the TF-gene regulation.
    S2
        Pandas dataframe containing network 2. Columsn are TFs and rows are
        regulated genes. Values are the weights of the TF-gene regulation.
    conditioning_val
        Soft thresholding value used to non TF-gene putative interactions

    """
    Y1 = _exp_df(c1)
    Y2 = _exp_df(c2)

    S1, conditioning_val = _network_df(n1, Y1, conditioning,
                                       weightNormalize=weightNormalize)
    S2, _ = _network_df(n2, Y2, conditioning, weightNormalize=weightNormalize)

    if n1 != n2:
        genes = list(set(np.concatenate([S1.index, S2.index])))
        tfs = list(set(np.concatenate([S1.columns, S2.columns])))
        S1 = S1.ix[genes, tfs].fillna(conditioning_val)
        S2 = S2.ix[genes, tfs].fillna(conditioning_val)
        Y1 = Y1.ix[S1.index, :].fillna(0)
        Y2 = Y2.ix[S2.index, :].fillna(0)
    else:
        Y1 = Y1.ix[S1.index, :].fillna(0)
        Y2 = Y2.ix[S2.index, :].fillna(0)

    X1 = Y1.ix[S1.columns, :]
    X2 = Y2.ix[S2.columns, :]
    _eval_indices(Y1, Y2, S1, S2)
    if null:
        Y1, Y2, X1, X2 = _shuffled_inputs(Y1, Y2, X1, X2, 
                seed=seed, shuffleGenes=shuffleGenes)
    return (Y1, Y2, X1, X2, S1, S2, conditioning_val)


def _gcp(x, r):
    """To calculate the graph contrained penalty term.

    Function that tensorflow uses to fold through each column of the matrix. It
    goes through each TF column vector, finds target genes, creates weights
    vector of the target genes, calculates pairwise difference of the weights,
    performs tanh tranformation of the differences, and then sums the
    difference.

    Parameters
    ----------
    x
        is the previous value
    r
        Weights of the TF matrix

    Returns
    -------
    float
        Value with sum of the weight difference between the target genes

    """
    ind = tf.where(tf.abs(r) > 0)
    vec = tf.expand_dims(tf.gather_nd(r, ind), 0)
    pc = tf.tanh(tf.transpose(vec)-vec)
    val = tf.divide(tf.reduce_sum(tf.abs(pc)), 2)
    return x+val


def run_model(Y, X, S, step, itr, log_itr, seed,
              l_reg=1e-4, g_reg=1e-4, stopthreshold=0.01, val=0,
              model='epee-gcl'):
    """To run sparse linear model.

    There are two sparse linear model implemented: lasso and
    graph-constrained-lasso. By default, method runs graph-constrained-lasso.

    Parameters
    ----------
    Y
        Pandas dataframe containing expression data. Rows are genes and
        columns are samples. Values are log(RPKM/TPM/FPKM + 1)
    X
        Pandas dataframe containing expression data of TFs. Rows are TFs and
        columns are samples. Values are log(RPKM/TPM/FPKM + 1)
    S
        Pandas dataframe containing network. Rows are genes and columns are
        TFs. Values are weight corresponding to the TF regulating a gene
    step
        learning rate for the optimizer
    log_itr
        Iterations to log the loss and percent change
    seed
        Setting the tensforflow random seed
    l_reg
        lasso regularization constant
    g_reg
        graph constrained regularization constant
    stopthreshold
        threshold when to stop learning the model if the loss change is 0.1
        between the previous iteration and the current interation
    val
        weight given to the not known TF-gene pairs
    model
        model to use for regulator and perturbed gene inference score

    Returns
    -------
    curr_y
        Inferred Y
    curr_w
        Inferred W
    loss_arr
        Loss per each iteration

    """
    tf.set_random_seed(seed)
    genes, samples = Y.shape
    regulators = X.shape[0]
    S_h = np.copy(S)
    S_h = np.float32(S_h > val)

    with tf.Graph().as_default():
        w = tf.Variable(
            tf.random_gamma([genes, regulators], alpha=1, beta=30,
                            dtype=tf.float32, seed=seed))

        if model == 'no-penalty':
            # least squares loss along with L1 regularization
            y = tf.matmul(tf.multiply(w, S), X)

            # loss
            loss = tf.reduce_mean(tf.square(Y - y))

        if model == 'epee-l':
            # least squares loss along with L1 regularization
            y = tf.matmul(tf.multiply(w, S), X)

            # loss
            loss = tf.reduce_mean(
                tf.square(Y - y))+tf.multiply(tf.reduce_sum(tf.abs(w)), l_reg)

        if model == 'epee-gcl':
            y = tf.matmul(tf.multiply(w, S), X)
            wa = tf.transpose(
                tf.multiply(tf.expand_dims(w, 1), tf.expand_dims(S_h, 1)))

            gc = tf.foldl(
                _gcp, wa, initializer=tf.constant(0, dtype=tf.float32),
                parallel_iterations=10, back_prop=False, swap_memory=True)

            loss = tf.reduce_mean(
                tf.square(Y-y)) + tf.multiply(
                    tf.reduce_sum(tf.abs(w)), l_reg) + tf.multiply(gc, g_reg)

        # optimizer
        optimizer = tf.train.AdamOptimizer(learning_rate=step, epsilon=1e-8)

        train = optimizer.minimize(loss)

        # training loop
        init = tf.global_variables_initializer()  # before starting init var
        sess = tf.Session()  # launch the graph
        sess.run(init)  # reset values to wrong

        loss_arr = []
        # outarr = []
        for s in range(itr):
            sess.run(train)
            if s % 100 == 0:
                curr_loss = sess.run(loss)
                if np.isnan(curr_loss):
                    raise RuntimeError('NAN value is computed for the loss.\
                    Make sure that the inputs are RPKM, TPM, FPKM normalized\
                    without any transformation (ie. log). Also can try to \
                    lower the learning rate.')
                loss_arr.append(curr_loss)
                if len(loss_arr) >= 2:
                    delta = loss_arr[-2] - loss_arr[-1]

                    logging.debug((s, curr_loss, delta))
                    if delta < 1:
                        logging.debug("TRANING FINISHED")
                        break
                else:
                    logging.debug((s, curr_loss))

        curr_w, curr_y, curr_loss = sess.run([w, y, loss])
    return (curr_y, curr_w, loss_arr)


def get_perturb_scores(Y1, y1, X1, w1, Y2, y2, X2, w2, S1, S2):
    """To get perturb scores.

    The funtion calculates the log likelihood ratio of error by swapping the
    weights between condition and error within same condition. The function
    returns sorted dataframe containing perturbation score.

    Parameters
    ----------
    Y1
        Pandas dataframe containing condition 1 expression data
    y1
        Pandas dataframe containing inferred condition 1 expression data
    X1
        Pandas dataframe containing TF expression across condition 1 samples
    w1
        Pandas dataframe containing inferred TF-gene weights for condition 1
    Y2
        Pandas dataframe containing condition 2 expression data
    y2
        Pandas dataframe containing inferred condition 2 expression data
    X2
        Pandas dataframe containing TF expression across condition 2 samples
    w2
        Pandas dataframe containing inferred TF-gene weights for condition 2
    S1
        Pandas dataframe containing network 1. Columns are TFs and rows are
        regulated genes. Values are the weights of the TF-gene regulation.
    S2
        Pandas dataframe containing network 2. Columsn are TFs and rows are
        regulated genes. Values are the weights of the TF-gene regulation.

    Returns
    -------
    scores_df_sorted
        Pandas dataframe containing perturb scores

    """
    err11 = np.square(Y1-y1).sum(axis=1)
    err22 = np.square(Y2-y2).sum(axis=1)
    err12 = np.square(Y1-np.dot(np.multiply(w2, S2), X1)).sum(axis=1)
    err21 = np.square(Y2-np.dot(np.multiply(w1, S1), X2)).sum(axis=1)

    base = err11+err22
    err = err21+err12
    scores = err/base

    # sort_scores_idx = sorted(range(len(scores)), key=lambda k: scores[k])
    scores_df = pd.DataFrame({'gene': scores.index, 'score': scores})
    scores_df_sorted = scores_df.sort_values(by='score', ascending=False)
    return scores_df_sorted


def get_summary_scoresdf(df, metric='sum'):
    """To calculate scores from multiple models.

    Each independent model generates perturb and regulatory score. The function
    calculates summary score for each perturb gene and assigns that score to
    the gene. The function returns dataframe with summary scores sorted.

    Parameters
    ----------
    df
        Pandas dataframe containing scores from multiple models
    metric
        Metric used to summarize the score from multiple models. 
        'sum', 'mean' and 'median' are valid options. Default = 'sum'

    Returns
    -------
    out_df_sort
        Pandas dataframe containing the summarized scores

    """
    if metric == 'median':
        df_score = df.iloc[:, 1:].median(axis=1)
    if metric == 'sum':
        df_score = df.iloc[:, 1:].sum(axis=1)
    if metric == 'mean':
        df_score = df.iloc[:, 1:].mean(axis=1)
    out_df = pd.DataFrame({'gene': df['gene'], 'score': df_score})
    out_df_sort = out_df.sort_values(by='score', ascending=False)
    out_df_sort.reset_index(inplace=True, drop=True)
    return out_df_sort


def get_weights_df(w, genes, tfs):
    """To name the rows and columns of the weight numpy ndarray object.

    Functions converts the numpy ndarray object to Pandas dataframe with
    labeled rows and columns.

    Parameters
    ----------
    w
        numpy ndarray containing the inferred TF-gene weights
    genes
        rownames of the w
    tfs
        column names of the w
    Returns
    -------
    df
        Pandas dataframe containing TF-gene inferred weights W

    """
    df = pd.DataFrame(w)
    df.columns = tfs
    df.index = genes
    return df


def get_regulator_scores(perturb_genes, W):
    """To provide regulator score given weights and perturb genes.

    Function calculates the regulator score given list of perturb genes and
    inferred W.

    Parameters
    ----------
    perturb_genes
        list of perturb genes
    W
        Pandas dataframe containing TF-gene inferred weights W

    Returns
    -------
    sort_scores_df
        Pandas dataframe containing the regulatory scores

    """
    # perturb_genes
    W_perturb = W.iloc[W.index.isin(perturb_genes), ]
    # not perturb_genes
    W_notperturb = W.iloc[~W.index.isin(perturb_genes), ]
    reg_score_outarr = []
    for reg in W.columns:
        num_score = np.sum(W_perturb[reg])
        not_reg_df = W_perturb.iloc[:, ~W_perturb.columns.isin([reg])]
        dem_score1 = not_reg_df.sum().sum()
        dem_score2 = np.sum(W_notperturb[reg])
        if dem_score1 == 0:
            dem_score1 = 1e-6
        if num_score == 0:
            score = 0
        else:
            score = (num_score/np.sqrt(dem_score1*dem_score2))
        reg_score_outarr.append(score)

    score_df = pd.DataFrame({'id': W.columns, 'score': reg_score_outarr})
    sort_scores_df = score_df.sort_values(by='score', ascending=False)
    sort_scores_df.reset_index(inplace=True, drop=True)
    return sort_scores_df


def get_diff_regulatory_activity(perturbed_genes, w1, w2, top_regs=10):
    """To provide differential regulator scores between condition.

    Function calculates the differential regulator score given list of perturb
    genes, and inferred W from the two conditions.

    Parameters
    ----------
    perturb_genes
        list of perturb genes
    w1
        Pandas dataframe containing TF-gene inferred weights W from condition 1
    w2
        Pandas dataframe containing TF-gene inferred weights W from condition 2
    top_regs
        Top number of differential regulators to select

    Returns
    -------
    reg_score
        Pandas dataframe containing the differential regulator scores
    diff_reg
        Pandas dataframe containing the top number of differential regulators

    """
    reg_score_w1 = get_regulator_scores(perturbed_genes, np.abs(w1))
    reg_score_w2 = get_regulator_scores(perturbed_genes, np.abs(w2))
    merge_reg_score = pd.merge(reg_score_w1, reg_score_w2, on='id')
    merge_reg_score.columns = ['id', 'w1', 'w2']
    merge_reg_score['w2-1'] = merge_reg_score['w2']-merge_reg_score['w1']
    diff_sorted_reg_score = merge_reg_score.sort_values(
        by='w2-1', ascending=False)
    diff_sorted_reg_score.reset_index(drop=True, inplace=True)
    reg_score = diff_sorted_reg_score[['id', 'w2-1']]
    diff_reg = pd.concat([reg_score.head(top_regs), reg_score.tail(top_regs)])
    diff_reg.reset_index(drop=True, inplace=True)
    return(reg_score, diff_reg)


def get_significant_scores(df, nullvalues, two_sided=False):
    """To calculate the score significance from emperical null CDF.

    Function generates emperical null CDF from the null values and calculate
    significance of the scores from the emperical CDF. Pvalues are corrected
    for the multiple hypothesis testing through Benjamini-Hochberg procedure.

    Parameters
    ----------
    df
        Pandas dataframe containing the perturb or regulator score
    nullvalues
        list of null values
    two_sided
        Boolean. Whether to use two-sided or one-sided test. Default = False.

    Returns
    -------
    df
        Pandas dataframe containing the pvalues

    """
    emcdf = sm.distributions.ECDF(nullvalues)
    median = np.median(nullvalues)
    pvals = []
    for score in df['score']:
        if two_sided:
            if score > median:
                pval = 1-emcdf(score)
            else:
                pval = emcdf(score)
        else:
            pval = 1-emcdf(score)
        pvals.append(pval)
    df['pvals'] = pvals
    df['fdr_bh'] = statsmodels.stats.multitest.multipletests(
        pvals, alpha=0.25, method='fdr_bh', returnsorted=False)[1]
    return df


def get_null_scores(dir):
    """To get null perturb scores.

    The functions generates null perturb score list.

    Parameters
    ----------
    dir
        Directory containing the perturbation scores from shuffled labels

    Returns
    -------
    null_scores_list
        list of null perturb scores

    """
    files = glob.glob('{}/pert*'.format(dir))
    null_scores = pd.DataFrame()
    for file in files:
        df = pd.read_csv(file, sep='\t', index_col=0)
        if null_scores.shape[0] == 0:
            null_scores = df
        else:
            null_scores = pd.merge(null_scores, df, on='gene')
    null_scores_list = null_scores.set_index('gene').values.flatten()
    return null_scores_list


def get_null_regscores(dir, perturbed_genes):
    """To get null regulatory scores.

    The function generates null regulatory scores from directory containing the
    inferred W's from the shuffled labels.

    Parameters
    ----------
    dir
        Directory containing the shuffled labels inferred W's.
    perturbed_genes
        list of perturbed genes used from the comparisions with true labels

    Returns
    -------
    null_score_list
        list of null regulator scores

    """
    files = glob.glob('{}/WS1*'.format(dir))
    null_regscore_df = pd.DataFrame()
    for file in files:
        seed = file.split('_')[-1].split('.')[0]
        ws1 = pd.read_csv(file, sep='\t', index_col=0)
        ws2 = pd.read_csv(file.replace('WS1', 'WS2'), sep='\t', index_col=0)
        null_regscore, _ = get_diff_regulatory_activity(
                                     perturbed_genes,
                                     ws1, ws2, top_regs=20)
        null_regscore.columns = ['gene', seed]
        if null_regscore_df.shape[0] == 0:
            null_regscore_df = null_regscore
        else:
            null_regscore_df = pd.merge(null_regscore_df, null_regscore,
                                        on='gene')
    null_score_list = null_regscore_df.set_index('gene').values.flatten()
    return null_score_list


def get_median_rank_position(actual, predicted, n, model):
    """To get median relative rank poistion.

    The method is used when ground truth reference is known and with query list
    want to known median rank position of the reference objects.

    Parameters
    ----------
    actual
        array of genes that are ground truth
    predicted
        array of genes that want to query
    n
        length of the query list
    model
        model used to get the query list

    Returns
    -------
    out
        median rank

    """
    predicted = pd.DataFrame(predicted).reset_index().set_index(0)
    predicted['index'] = (predicted['index']/n)*100
    rank = []
    for r in predicted.ix[actual, :].iterrows():
        val = r[1][0]
        rank.append(val)
    out = (np.median(rank), model)
    return out
