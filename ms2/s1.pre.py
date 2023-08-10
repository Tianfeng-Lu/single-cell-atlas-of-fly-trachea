import os,sys
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import pearsonr


df = pd.read_csv('all.data.csv', index_col = 0)


G = ['SB','DT','DB','ASP','TC','LT','VB','PC','GB']
group = [['PC', 'PC_S1'], ['SB', 'SB_M1'], ['ASP', 'ASP_S1'], ['ASP', 'ASP_S2'], ['ASP', 'ASP_S3']]

dfo = pd.DataFrame()
for i in group:
    gene = []
    other = [x for x in G if x not in i]
    if i[0] != 'ASP':
        df1 = df[i]
    else:
        df1 = df[[i[0]]]
    df2 = df[other]
    for j in df1.index:
        v2 = list(df2.loc[j])
        if i[0] != 'ASP':
            v1 = list(df1.loc[j])
            if np.mean(v1) < np.mean(v2) * 2:
                continue
            p = stats.ttest_ind(v1, v2)
            if p.__getattribute__("pvalue") > 0.05:
                continue
        else:
            if list(df1.loc[j])[0] < np.median(v2) * 2:
                continue
        gene.append(j)
    dfs = df.loc[gene]
    value = []
    for j in G:
        value.append(pearsonr(list(dfs[i[1]]), list(dfs[j]))[0])
    dfo[i[1]] = value

dfo.index = G
dfo.to_csv('all.pearson.csv')

