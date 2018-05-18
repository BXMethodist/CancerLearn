import os, pandas as pd, numpy as np
from sklearn import preprocessing
from scipy import stats


def center_normalization(table):
    return preprocessing.StandardScaler().fit(table)

def preprocessing_table(scaler, table):
    new_table_data = scaler.transform(table)
    new_table = pd.DataFrame(new_table_data, index=table.index, columns=table.columns)
    return new_table

def label_label(table, label_table):
    result = [0]*table.shape[0]
    table['label'] = result
    print label_table['label'].unique()
    for i in label_table['label'].unique():
        cur_label = label_table[label_table['label']==i]
        table.ix[table.index.isin(cur_label.index), 'label'] = i
    return table

def logP_wilcoxon(groupA, groupB, bags=1):
    """
    return the negative log P value for two groups
    :param groupA:
    :param groupB:
    :return:
    """
    p_values = []

    for i in range(bags):
        # cur_groupA = np.random.choice(groupA, int(len(groupA)*0.75), replace=True)
        # cur_groupB = np.random.choice(groupB, int(len(groupB)*0.75), replace=True)
        cur_groupA = groupA
        cur_groupB = groupB
        try:
            rank_diff, p1 = stats.mannwhitneyu(cur_groupA, cur_groupB, alternative='less')
            rank_diff, p2 = stats.mannwhitneyu(cur_groupA, cur_groupB, alternative='less')
            if p1 > p2:
                p = p2
            else:
                p = p1
            p_values.append(np.log10(p))
        except:
            p_values.append(0)

    return np.mean(p_values)


def tau(final_df):
    cdf = final_df + abs(final_df.min().min())
    return (1 - cdf.divide(cdf.max(axis=1), axis=0)).sum(axis=1) / (cdf.shape[1] - 1)


def quantileNormalize(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df

def submit_pbs(cmds, outputname, queue='highmem', time='32:00:00', ppn='1', mem="16000mb"):
    """
    :param cmd: command line
    :param outputname: output folder name
    :return:
    """
    pbs = open(outputname + ".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N " + outputname + '\n')
    pbs.write("#PBS -q "+queue+"\n")
    pbs.write("#PBS -m e\n")
    pbs.write("#PBS -M bxia@houstonmethodist.org\n")
    pbs.write("#PBS -l walltime="+time+"\n")
    pbs.write("#PBS -l nodes=1:ppn="+ppn+"\n")
    pbs.write("#PBS -l pmem="+mem+"\n")
    pbs.write("cd " + os.getcwd() + "\n")
    pbs.write("module load python/2.7.11\n")
    pbs.write("module load R/3.2.1\n")
    for cmd in cmds:
        pbs.write(cmd + "\n")
    pbs.close()
    os.system('qsub ' + outputname + ".pbs")