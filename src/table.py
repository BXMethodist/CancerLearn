"""
This module take a table with best parameter combinations, and generate the table with features for each genes.
"""
import os, pandas as pd, numpy as np, random
from collections import defaultdict
from copy import deepcopy
from multiprocessing import Process, Queue
from itertools import permutations
from scipy.stats import norm

def random_control_genes(CIG_gene_df, all_genes, random_seed, number_genes=None):
    """
    random select genes for each cell type

    :param CIG_gene_table: table path
    :param all_genes: all gene GTF
    :param random_seed: random seed number to make sure each return the same random gene set
    :return:
    """

    results = []
    random.seed(random_seed)
    cur_number = CIG_gene_df.shape[0]
    cur_candidates = set()
    for gene in all_genes:
        if gene not in CIG_gene_df['gene'].unique():
            cur_candidates.add(gene)

    # print 'ISOC1' in cur_candidates, 'why CIG in?'

    if number_genes is not None:
        cur_negative = random.sample(cur_candidates, number_genes)
    else:
        cur_negative = random.sample(cur_candidates, cur_number)
    for i in range(len(cur_negative)):
        results.append([cur_negative[i], 'CT', 0])

    negative_control_genes_df = pd.DataFrame(results)
    negative_control_genes_df.columns = CIG_gene_df.columns
    return negative_control_genes_df

def df_to_index_danpos(df, bin=3000):
    results = {}
    f = open(df, 'r')
    for line in f.readlines()[1:]:
        line = line.split()
        if len(line) <=1:
            continue
        t = (line[0],
             int(line[1]),
             int(line[2]),
             int(line[4]),
             float(line[5]),
             float(line[6]),
             )

        if t[0] not in results:
            results[t[0]] = {}

        for i in range(t[1]/bin, t[2]/bin+1):
            if i in results[t[0]]:
                results[t[0]][i].add(t)
            else:
                results[t[0]][i] = set()
                results[t[0]][i].add(t)
    f.close()
    return results

def df_to_index_sk(df, bin=3000):
    results = {}
    # print df
    f = open(df, 'r')

    for line in f.readlines()[1:]:
        line = line.strip().split(',')
        t = (line[0],
             int(line[1]),
             int(line[2]),
             int(line[4]),
             float(line[5]),
             float(line[6]),
             float(line[8]),
             float(line[9]),
             )

        if t[0] not in results:
            results[t[0]] = {}

        for i in range(t[1]/bin, t[2]/bin+1):
            if i in results[t[0]]:
                results[t[0]][i].add(t)
            else:
                results[t[0]][i] = set()
                results[t[0]][i].add(t)
    f.close()
    return results

def get_range_absolute(gene_list, all_gene_GTF, left_distance, right_distance, TSS_pos, TTS_pos):
    """

    :param all_gene_GTF:
    :param gene_list:
    :param left_distance: upstream is -, downstream is +
    :param right_distance: upstream is -, downstream is +
    :param left_pos:
    :param right_pos:
    :return:
    """

    if gene_list is not None:
        cur_df = all_gene_GTF[all_gene_GTF['hg19.kgXref.geneSymbol'].isin(gene_list)]
    else:
        cur_df = all_gene_GTF

    positive_df = cur_df[cur_df['hg19.knownGene.strand'] == '+'].copy()
    negative_df = cur_df[cur_df['hg19.knownGene.strand'] == '-'].copy()

    if TSS_pos == 'TSS' and TTS_pos == 'TSS':
        positive_df['left_range'] = positive_df['hg19.knownGene.txStart'] + left_distance
        positive_df.loc[positive_df['left_range'] < 0, 'left_range'] = 0
        positive_df['right_range'] = positive_df['hg19.knownGene.txStart'] + right_distance

        negative_df['right_range'] = negative_df['hg19.knownGene.txEnd'] - left_distance
        negative_df['left_range'] = negative_df['hg19.knownGene.txEnd'] - right_distance
        negative_df.loc[negative_df['left_range'] < 0, 'left_range'] = 0

    elif TSS_pos == 'TSS' and TTS_pos == 'TTS':
        positive_df['left_range'] = positive_df['hg19.knownGene.txStart'] + left_distance
        positive_df.loc[positive_df['left_range'] < 0, 'left_range'] = 0
        positive_df['right_range'] = positive_df['hg19.knownGene.txEnd'] + right_distance

        negative_df['right_range'] = negative_df['hg19.knownGene.txEnd'] - left_distance
        negative_df['left_range'] = negative_df['hg19.knownGene.txStart'] - right_distance
        negative_df.loc[negative_df['left_range'] < 0, 'left_range'] = 0

    elif TSS_pos == 'MID':
        positive_df['MID'] = (positive_df['hg19.knownGene.txStart'] + positive_df['hg19.knownGene.txEnd'])/2
        negative_df['MID'] = (negative_df['hg19.knownGene.txStart'] + negative_df['hg19.knownGene.txEnd'])/2

        positive_df['left_range'] = positive_df['MID'] + left_distance
        positive_df[positive_df['left_range'] < 0] = 0
        positive_df['right_range'] = positive_df['MID'] + right_distance

        negative_df['right_range'] = negative_df['MID'] - left_distance
        negative_df['left_range'] = negative_df['MID'] - right_distance
        negative_df[negative_df['left_range'] < 0] = 0
    positive_df['length'] = positive_df['hg19.knownGene.txEnd'] - positive_df['hg19.knownGene.txStart']
    negative_df['length'] = negative_df['hg19.knownGene.txEnd'] - negative_df['hg19.knownGene.txStart']

    new_df = positive_df.append(negative_df)
    # print new_df.shape, cur_df.shape, len(new_df['hg19.kgXref.geneSymbol'].unique())

    if len(new_df['hg19.kgXref.geneSymbol'].unique()) != len(all_gene_GTF['hg19.kgXref.geneSymbol'].unique()):
        new_df.to_csv('wrong_range.csv')
        # print new_df[new_df['hg19.kgXref.geneSymbol'] ==0]
        # print all_gene_GTF[all_gene_GTF['hg19.kgXref.geneSymbol'] == 0]

    result_df = new_df[['hg19.kgXref.geneSymbol', 'hg19.knownGene.chrom', 'left_range', 'right_range', 'length']]
    result_df.columns = ['gene', 'chr', 'left_range', 'right_range', 'length']

    return result_df

def get_stats(gene_df, df_path, criteria, bin=3000, df_function=df_to_index_danpos):
    """
    This function will select the target gene's peaks stats from the dataframe
    :param gene_df: gene:range dictionary, get the best result for each gene from transcripts range
    :param df: dataframe contain peaks from certain cutoff
    :return:
    """
    # print 'get stats', len(gene_df['gene'].unique())
    if criteria != 'skewness' and criteria != 'kurtosis':
        table_dict = df_function(df_path)
    else:
        df_function = df_to_index_sk
        table_dict = df_function(df_path)

    results = defaultdict(float)

    for k in range(gene_df.shape[0]):
        gene_name = gene_df.iloc[k, 0]

        chr_name, start, end, length = gene_df.iloc[k, 1], gene_df.iloc[k, 2], gene_df.iloc[k, 3], gene_df.iloc[k, 4]
        ## Here is the problem, danpos selector will consider the entire overlapped peaks
        ## The other approach is using self designed peak calling, to make sure each parameter will return different value
        cur_table = set()

        if end < start:
            mid = (start + end) / 2
            start = mid
            end = mid

        for i in range(int(start/bin), int(end/bin) + 1):
            if chr_name in table_dict and i in table_dict[chr_name]:
                table = table_dict[chr_name][i]
                cur_table = cur_table.union(table)

        if len(cur_table) == 0:
            continue

        selected_table = []
        for t in cur_table:
            if start < t[1] < end:
                selected_table.append(t)
            elif start < t[2] < end:
                selected_table.append(t)
            elif t[1] <= start and end <= t[2]:
                selected_table.append(t)

        if len(selected_table) == 0:
            continue

        cur_df = pd.DataFrame(list(selected_table))

        if cur_df.shape[1] == 6:
            cur_df.columns = ['chr',
                          'start',
                          'end',
                          'width_above_cutoff',
                          'total_signal',
                          'height',]
        else:
            cur_df.columns = ['chr',
                              'start',
                              'end',
                              'width_above_cutoff',
                              'total_signal',
                              'height',
                              'skewness',
                              'kurtosis']

        if criteria == 'total_width':
            cur_col = cur_df['end'] - cur_df['start']
            cur_value = cur_col.sum()
        elif criteria == 'height':
            cur_value = cur_df['height'].max()
        elif criteria == 'single_width':
            cur_col = cur_df['end'] - cur_df['start']
            cur_value = cur_col.max()
        elif criteria == 'total_signal':
            cur_value = cur_df['total_signal'].sum()
        elif criteria == 'single_signal':
            cur_value = cur_df['total_signal'].max()
        elif criteria == 'coverage':
            cur_value = (cur_df['end'] - cur_df['start']).sum()*1.0/length


        #
        # # This is for kurtosis and skewness
        elif cur_df.shape[0] > 0 and criteria == 'skewness' and 'skewness' in cur_df.columns:
            cur_value = cur_df.ix[cur_df['total_signal'].argmax(),'skewness']
        elif cur_df.shape[0] > 0 and criteria == 'kurtosis' and 'kurtosis' in cur_df.columns:
            cur_value = cur_df.ix[cur_df['total_signal'].argmax(), 'kurtosis']


        if cur_value > results[gene_name] and criteria != 'skewness' and criteria != 'kurtosis':
            results[gene_name] = cur_value
        # this is for kurtosis and skewness

        elif criteria == 'kurtosis':
            if abs(cur_value) > abs(results[gene_name]):
                results[gene_name] = cur_value
        elif criteria == 'skewness':
            if abs(cur_value) > results[gene_name]:
                results[gene_name] = abs(cur_value)

    final = []

    for gene_name in gene_df['gene'].unique():
        final.append((gene_name, results[gene_name]))
    print len(final)
    return final

def CIG_selecter(CIG_gene_df, non_CIG_gene_df, all_gene_GTF, up_stream_distance, down_stream_distance, all_dfs, cutoff, marker, criteria,
                 TSS_pos, TTS_pos):
    """
    get the genes status and return a data frame with two columns, gene name and criteria.
    :param CIG_gene_df:
    :param non_CIG_gene_df:
    :param all_gene_GTF:
    :param up_stream_distance:
    :param widow_size:
    :param all_dfs: dictionary of dictionary, cell type and cutoff
    :param cutoff:
    :return:
    """
    CIG_gene_list = list(CIG_gene_df['gene'].values)
    non_CIG_gene_list = list(non_CIG_gene_df['gene'].values)



    CIG_gene_ranges = get_range_absolute(CIG_gene_list, all_gene_GTF, up_stream_distance, down_stream_distance,
                                         TSS_pos,
                                         TTS_pos)
    non_CIG_gene_ranges = get_range_absolute(non_CIG_gene_list, all_gene_GTF, up_stream_distance, down_stream_distance,
                                             TSS_pos,
                                             TTS_pos)

    CIG_results = []
    non_CIG_results = []

    cur_CIG_gene_list = list(CIG_gene_df['gene'].values)
    cur_CIG_gene_range = CIG_gene_ranges[CIG_gene_ranges['gene'].isin(cur_CIG_gene_list)]

    cur_non_CIG_gene_list = list(non_CIG_gene_df['gene'].values)
    cur_non_CIG_gene_range = non_CIG_gene_ranges[non_CIG_gene_ranges['gene'].isin(cur_non_CIG_gene_list)]

    cur_df = all_dfs[cutoff]

    cur_CIG_result = get_stats(cur_CIG_gene_range, cur_df, criteria)
    cur_non_CIG_result = get_stats(cur_non_CIG_gene_range, cur_df, criteria)

    CIG_results += cur_CIG_result
    non_CIG_results += cur_non_CIG_result


    CIG_results_df = pd.DataFrame(CIG_results)

    if TTS_pos == 'TTS':
        CIG_results_df.columns = ['gene', marker + '_' + criteria + '_' + 'genebody']
    else:
        CIG_results_df.columns = ['gene', marker + '_' + criteria]
    CIG_results_df=CIG_results_df.set_index(['gene'])

    non_CIG_results_df = pd.DataFrame(non_CIG_results)
    if TTS_pos == 'TTS':
        non_CIG_results_df.columns = ['gene', marker + '_' + criteria + '_' + 'genebody']
    else:
        non_CIG_results_df.columns = ['gene', marker +'_'+ criteria]
    non_CIG_results_df = non_CIG_results_df.set_index(['gene'])

    return CIG_results_df, non_CIG_results_df

def CIG_selecter_all(all_gene_GTF, up_stream_distance, down_stream_distance, all_dfs, cutoff, criteria, cell_type,
                     TSS_pos, TTS_pos):
    """
    get the genes status and return a data frame with two columns, gene name and criteria.
    :param CIG_gene_df:
    :param non_CIG_gene_df:
    :param all_gene_GTF:
    :param up_stream_distance:
    :param widow_size:
    :param all_dfs: dictionary of dictionary, cell type and cutoff
    :param cutoff:
    :return:
    """
    all_gene_ranges = get_range_absolute(None, all_gene_GTF, up_stream_distance, down_stream_distance,
                                         TSS_pos, TTS_pos)

    cur_df = all_dfs[cell_type][cutoff]

    # overlap selecter
    all_gene_results = get_stats(all_gene_ranges, cur_df, criteria, cell_type)

    all_gene_results_df = pd.DataFrame(all_gene_results)

    all_gene_results_df.columns = ['gene', criteria, 'cell_type']
    all_gene_results_df = all_gene_results_df[['gene', criteria]]
    return all_gene_results_df

def get_training_table(all_gene_GTF_path='/archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.xls',
                       CIG_gene_df_path='/archive/tmhbxx3/CancerLearn/genelist/final_cancer_genes.xlsx',
                       control_gene_df_path='/archive/tmhbxx3/CancerLearn/genelist/Control_genes.xlsx',
                       best_parameter_path='../ref_data/parameters_default.csv',
                       process=8, size=500, random_seed=1,
                       out_name=None):
    ## get the CIG genesets and nonCIG genesets
    all_gene_GTF = pd.read_csv(all_gene_GTF_path, sep='\t', dtype={'hg19.kgXref.geneSymbol': str})
    all_gene_GTF['hg19.kgXref.geneSymbol'] = all_gene_GTF['hg19.kgXref.geneSymbol'].str.upper()
    all_gene_GTF = all_gene_GTF[~(all_gene_GTF['hg19.kgXref.geneSymbol'].str.contains('METAZOA_SRP'))]
    all_gene_GTF = all_gene_GTF[~(all_gene_GTF['hg19.kgXref.geneSymbol'].str.contains('D28408'))]
    all_gene_GTF = all_gene_GTF[~(all_gene_GTF['hg19.kgXref.geneSymbol'].str.contains('MIR'))]
    all_gene_GTF = all_gene_GTF[~(all_gene_GTF['hg19.kgXref.geneSymbol'].str.contains('-AS'))]
    all_gene_GTF = all_gene_GTF[~(all_gene_GTF['hg19.kgXref.geneSymbol'].str.contains('LINC'))]
    all_gene_GTF = all_gene_GTF[~(all_gene_GTF['hg19.kgXref.geneSymbol'].str.startswith('LOC'))]
    all_gene_GTF = all_gene_GTF[~(all_gene_GTF['hg19.kgXref.geneSymbol'].str.startswith('TRNA_'))]
    all_gene_GTF = all_gene_GTF[~(all_gene_GTF['hg19.kgXref.geneSymbol'].str.contains(r"[J][A][0-9]+"))]
    all_gene_GTF = all_gene_GTF[~(all_gene_GTF['hg19.kgXref.geneSymbol'].str.contains(r"[D][Q][0-9]+"))]
    all_gene_GTF = all_gene_GTF[~(all_gene_GTF['hg19.kgXref.geneSymbol'].str.contains(r"[J][B][0-9]+"))]
    all_gene_GTF = all_gene_GTF[~(all_gene_GTF['hg19.kgXref.geneSymbol'].str.contains(r"[H][V][0-9]+"))]


    control_gene_df = pd.read_excel(control_gene_df_path)
    control_gene_df['gene'] = control_gene_df['gene'].str.upper()
    control_genes = list(control_gene_df['gene'].unique())

    if CIG_gene_df_path.endswith('.csv'):
        CIG_gene_df = pd.read_csv(CIG_gene_df_path, dtype={'gene': str})
    elif CIG_gene_df_path.endswith('.xlsx'):
        CIG_gene_df = pd.read_excel(CIG_gene_df_path)
    CIG_gene_df['gene'] = CIG_gene_df['gene'].str.upper()

    ## get the full table
    all_genes = set(all_gene_GTF['hg19.kgXref.geneSymbol'].unique())
    cur_candidates = set()
    for gene in all_genes:
        if gene not in CIG_gene_df['gene'].unique():
            cur_candidates.add(gene)
    size = len(cur_candidates)
    # print CIG_gene_df.shape, 'CIG genes'

    ## for grid optimization
    # non_CIG_gene_df = random_control_genes(CIG_gene_df, exclude_list, all_genes, random_seed=1, number_genes=size, candidates=candidates)
    # non_CIG_gene_df['gene'] = non_CIG_gene_df['gene'].str.upper()
    # non_CIG_gene_df.to_csv('non_CIG_control_v6_grid.csv', index=None, encoding='utf-8')

    ## for lg
    non_CIG_gene_df = random_control_genes(CIG_gene_df, all_genes, random_seed=random_seed, number_genes=size)
    non_CIG_gene_df['gene'] = non_CIG_gene_df['gene'].str.upper()
    non_CIG_gene_df.to_csv('cancer_control'+str(size)+'_'+str(random_seed)+'.csv', index=None, encoding='utf-8')

    # return
    # non_CIG_gene_df = pd.read_csv('non_CIG_control.csv', dtype={'gene': str})
    # non_CIG_gene_df['gene'] = non_CIG_gene_df['gene'].str.upper()

    # get feature parameter table
    best_feature_df = pd.read_csv(best_parameter_path)

    # get all tables:
    all_tables = {}
    all_sk_tables = {}

    markers = os.listdir('/archive/tmhbxx3/CancerLearn/peaks/CD4/')

    for m in markers:
        if m not in all_tables:
            all_tables[m] = defaultdict(dict)
            dfs_path = '/archive/tmhbxx3/CancerLearn/peaks/CD4/'+m+'/'
            dfs = [x for x in os.listdir(dfs_path) if x.endswith('.xls')]

            for table_name in dfs:
                info = table_name.split('_')
                # print info
                cutoff = float(info[-1][:-4])
                all_tables[m][cutoff] = dfs_path+table_name

        if m not in all_sk_tables:
            all_sk_tables[m] = defaultdict(dict)
            dfs_path = '/archive/tmhbxx3/CancerLearn/peaks/CD4_sk/'+m+'/'
            dfs = [x for x in os.listdir(dfs_path) if x.endswith('.csv')]

            for table_name in dfs:
                info = table_name.split('_')
                cutoff = float(info[-1][:-4])
                all_sk_tables[m][cutoff] = dfs_path + table_name

    chunks = []
    cur_index = 0
    reminder = best_feature_df.shape[0] % process
    chunk_size = best_feature_df.shape[0] / process
    for i in range(process):
        if reminder > 0:
            chunks.append(best_feature_df[cur_index + i * chunk_size:cur_index + (i + 1) * chunk_size + 1])
            cur_index += 1
            reminder -= 1
        else:
            chunks.append(best_feature_df[cur_index + i * chunk_size: cur_index + (i + 1) * chunk_size])

    total_chunk_size = 0
    for chunk in chunks:
        total_chunk_size += len(chunk)
    if total_chunk_size != best_feature_df.shape[0]:
        print 'multiple processes chunk size is not correct'
        return None

    # print chunks
    queue = Queue()
    processes = []

    for i in range(process):
        cur_chunk = chunks[i]
        p = Process(target=training_table_process,
                    args=(queue, cur_chunk, CIG_gene_df, non_CIG_gene_df, all_gene_GTF, all_tables, all_sk_tables))
        processes.append(p)
        p.start()

    final_feature_tables = []

    for i in range(process):
        cur_final_feature_df = queue.get()
        final_feature_tables.append(cur_final_feature_df)
    final_feature_df = pd.concat(final_feature_tables, axis=1)

    for p in processes:
        p.join()

    final_feature_df.index.name = 'gene'
    if out_name is None:
        final_feature_df.to_excel('Cancer_train_500OGTSG.xlsx')
    else:
        final_feature_df.to_excel(out_name)


def training_table_process(queue, best_feature_df, CIG_gene_df, non_CIG_gene_df, all_gene_GTF, all_tables, all_sk_tables):
    cur_final_feature_df = pd.DataFrame(index=list(CIG_gene_df.gene.values) + list(non_CIG_gene_df.gene.values))
    # print best_feature_df.columns
    for i in range(best_feature_df.shape[0]):
        marker = best_feature_df.iloc[i, 0]
        feature = best_feature_df.iloc[i, 1]
        criteria = feature.replace('_genebody', '')

        start = best_feature_df.iloc[i, 2]
        end = best_feature_df.iloc[i, 3]
        height = best_feature_df.iloc[i, 4]
        print marker, feature, start, end, height
        if feature.find('genebody') != -1:
            TSS, TTS = 'TSS', 'TTS'
        else:
            TSS, TTS = 'TSS', 'TSS'

        if feature.find('kurtosis') != -1 or feature.find('skewness') != -1:
            option = True
        else:
            option = False

        if option:
            cur_stat_dfs = all_sk_tables[marker]
        else:
            cur_stat_dfs = all_tables[marker]

        cur_CIG_results_df, cur_non_CIG_results_df = CIG_selecter(CIG_gene_df, non_CIG_gene_df, all_gene_GTF,
                                                                  start, end,
                                                                  cur_stat_dfs, height, marker, criteria,
                                                                  TSS, TTS)
        cur_result_df = cur_CIG_results_df.append(cur_non_CIG_results_df)

        print len(cur_final_feature_df.index), len(cur_final_feature_df.index.unique())
        print len(cur_result_df.index), len(cur_result_df.index.unique())
        cur_final_feature_df[marker + '_' + feature] = cur_result_df[marker + '_' + feature]
    # cur_final_feature_df.to_csv(str(os.getpid())+'training_table.csv')
    queue.put(cur_final_feature_df)
    return


# def get_real_table(target_cell_type, all_gene_GTF_path='/archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.xls',
#                        best_parameter_path='best_parameters_CIG_6th_bo_1xControl.csv',
#                        process=8,
#                     out_name=None):
#     ## get the CIG genesets and nonCIG genesets
#     all_gene_GTF = pd.read_csv(all_gene_GTF_path, sep='\t', dtype={'hg19.kgXref.geneSymbol': str})
#
#     # all_gene_GTF = pd.read_csv('hg19.GREATgene2UCSCknownGenes.table.xls', sep='\t')
#     all_gene_GTF['hg19.kgXref.geneSymbol'] = all_gene_GTF['hg19.kgXref.geneSymbol'].str.upper()
#
#     all_genes = set(all_gene_GTF['hg19.kgXref.geneSymbol'].unique())
#
#     # get feature parameter table
#     best_feature_df = pd.read_csv(best_parameter_path)
#
#     # get all tables:
#     all_tables = {}
#     all_sk_tables = {}
#
#     # ['h3k4me3', 'h3k27me3', 'h3k4me1', 'h3k27ac']
#     for m in ['h3k4me3_qn', 'h3k27me3_qn', 'h3k4me1_qn', 'h3k27ac_qn']:
#         if m not in all_tables:
#             all_tables[m] = defaultdict(dict)
#             dfs_path = '/home/tmhbxx3/scratch/CIG/' + m + '_peaks/pooled/'
#             # dfs_path = '/home/tmhbxx3/scratch/CIG/test_peaks/' + m + '_peaks/pooled/'
#             dfs = [x for x in os.listdir(dfs_path) if x.endswith('.xls')]
#
#             for table_name in dfs:
#                 info = table_name.split('_')
#                 cell_type = info[2]
#                 cutoff = float(info[-1][:-4])
#                 # print cell_type, cutoff
#                 all_tables[m][cell_type][cutoff] = dfs_path+table_name
#
#         if m not in all_sk_tables:
#             all_sk_tables[m] = defaultdict(dict)
#             dfs_path = '/home/tmhbxx3/scratch/CIG/test/' + m + '_sk_peaks/'
#             dfs = [x for x in os.listdir(dfs_path) if x.endswith('.csv')]
#
#             for table_name in dfs:
#                 info = table_name.split('_')
#                 cell_type = info[0]
#                 cutoff = float(info[-1][:-4])
#                 all_sk_tables[m][cell_type][cutoff] = dfs_path + table_name
#
#     chunks = []
#     cur_index = 0
#     reminder = best_feature_df.shape[0] % process
#     chunk_size = best_feature_df.shape[0] / process
#     for i in range(process):
#         if reminder > 0:
#             chunks.append(best_feature_df[cur_index + i * chunk_size:cur_index + (i + 1) * chunk_size + 1])
#             cur_index += 1
#             reminder -= 1
#         else:
#             chunks.append(best_feature_df[cur_index + i * chunk_size: cur_index + (i + 1) * chunk_size])
#
#     total_chunk_size = 0
#     for chunk in chunks:
#         total_chunk_size += len(chunk)
#     if total_chunk_size != best_feature_df.shape[0]:
#         print 'multiple processes chunk size is not correct'
#         return None
#
#     queue = Queue()
#     processes = []
#
#     for i in range(process):
#         cur_chunk = chunks[i]
#         p = Process(target=real_table_process,
#                     args=(queue, target_cell_type, cur_chunk, all_gene_GTF, all_tables, all_sk_tables))
#         processes.append(p)
#         p.start()
#
#     final_feature_tables = []
#
#     for i in range(process):
#         cur_final_feature_df = queue.get()
#         final_feature_tables.append(cur_final_feature_df)
#     final_feature_df = pd.concat(final_feature_tables, axis=1)
#
#     for p in processes:
#         p.join()
#
#     index_cols = []
#     index_col = False
#     for ci in range(len(final_feature_df.columns)):
#         if final_feature_df.columns[ci].startswith('gene'):
#             if not index_col:
#                 index_cols.append(ci)
#                 index_col = True
#         else:
#             index_cols.append(ci)
#     final_feature_df = final_feature_df.iloc[:, index_cols]
#     columns = list(final_feature_df.columns)
#     columns[0] = 'gene_id'
#     final_feature_df.columns = columns
#
#     if out_name is None:
#         final_feature_df.to_csv(target_cell_type+'_' + best_parameter_path[:-4]+'_'+'real_table.csv', index=None)
#     else:
#         final_feature_df.to_csv(out_name, index=None)
#
# def real_table_process(queue, target_cell_type, best_feature_df, all_gene_GTF, all_tables, all_sk_tables):
#     cur_final_feature_df = None
#     # print best_feature_df.columns
#     for i in range(best_feature_df.shape[0]):
#         marker = best_feature_df.iloc[i, 0] + '_qn'
#         feature = best_feature_df.iloc[i, 1]
#         print feature
#         criteria = feature.replace('_genebody', '')
#
#         start = best_feature_df.iloc[i, 2]
#         end = best_feature_df.iloc[i, 3]
#         height = best_feature_df.iloc[i, 4]
#         # print marker, feature, start, end, height
#         if feature.find('genebody') != -1:
#             TSS, TTS = 'TSS', 'TTS'
#         else:
#             TSS, TTS = 'TSS', 'TSS'
#
#         if feature.find('kurtosis') != -1 or feature.find('skewness') != -1:
#             option = True
#         else:
#             option = False
#
#         if option:
#             cur_stat_dfs = all_sk_tables[marker]
#         else:
#             cur_stat_dfs = all_tables[marker]
#
#         all_gene_results_df = CIG_selecter_all(all_gene_GTF, start, end, cur_stat_dfs, height,
#                                                criteria, target_cell_type, TSS, TTS)
#
#         all_gene_results_df.columns = ['gene'+str(i), marker + '_' + feature]
#         if cur_final_feature_df is None:
#             cur_final_feature_df = all_gene_results_df.copy()
#         else:
#             cur_final_feature_df[marker + '_' + feature] = all_gene_results_df[all_gene_results_df.columns[1]]
#
#     # cur_final_feature_df.to_csv(str(os.getpid())+'_'+target_cell_type+'_'+'real_table.csv')
#     queue.put(cur_final_feature_df)
#     return

def add_RNA_expression(pathname, outname):
    gtf = pd.read_csv('hg19.ucscgenes.knowngene.xls', sep='\t', dtype={'hg19.kgXref.geneSymbol': str})

    d = 'genes_rankproduct.tsv'

    df = pd.read_csv(d, sep='\t', dtype={'gene_id': str})
    df = df.set_index(['gene_id'])
    df.index = df.index.str.upper()
    files = [x for x in os.listdir('.') if x.endswith(pathname)]
    print files
    for f in files[-1:]:
        target_df = pd.read_csv(f)
        columns = list(target_df.columns)
        columns[0] = 'gene_id'
        target_df.columns = columns
        results = []
        for i in range(target_df.shape[0]):
            if f.find('HUVEC') != -1:
                cell_type = 'HUVEC'
            else:
                cell_type = target_df.ix[i, 'cell_type']
            gene_id = target_df.ix[i, 'gene_id']
            try:
                value = df.ix[gene_id, cell_type+'_FPKM']
                if isinstance(value, float):
                    results.append(df.ix[gene_id, cell_type+'_FPKM'])
                else:
                    results.append(df.ix[gene_id, cell_type + '_FPKM'].max())
            except:
                print gene_id
                results.append(0)

        target_df['Rank_Product'] = results
        target_df.to_csv(outname, index=False)

def rank_product(table):
    # print table.product(axis=1)
    return np.power(table.product(axis=1), 1.0/table.shape[1])

def rank_row(table):
    result = table.rank(axis=0)
    result.fillna(result.shape[0], inplace=True)
    return result

def fold_change_rank(table, column):
    columns = [c for c in table.columns if c != column]
    cur_table = table[columns].copy()
    for c in cur_table.columns:
        cur_table[c] = cur_table[c]/table[column]
    cur_table.fillna(0, inplace=True)
    # print cur_table
    return rank_row(cur_table)

def rank_product_meta(table, iterations, number_permutated_samples=5):
    results = pd.DataFrame(index=table.index)
    for column in table.columns:
        print column
        cur_columns = [x for x in table.columns if x != column]
        cur_rank_table = fold_change_rank(table, column)
        # print cur_rank_table
        # print cur_rank_table
        # print cur_rank_table
        cur_rank_product = rank_product(cur_rank_table)
        # cur_rank = rank_row(cur_rank_product)
        cur_permutations = pd.DataFrame(index=table.index)

        cur_permutated_columns = permutations(list(cur_columns), number_permutated_samples)

        i = 0
        # print cur_rank_table
        for cur_p in cur_permutated_columns:
            cur_p_rank_table = cur_rank_table[list(cur_p)]
            cur_p_rank_product = rank_product(cur_p_rank_table)
            cur_p_rank_product[cur_p_rank_product <= cur_rank_product] = 1
            cur_p_rank_product[cur_p_rank_product > cur_rank_product] = 0
            cur_permutations[i]=cur_p_rank_product
            i += 1
            # print i
            if i >= iterations:
                break
        p_values = cur_permutations.mean(axis=1)

        p_values = 1 - norm.cdf(norm.ppf(1- p_values))

        results[column] = cur_rank_product
    results.to_csv('genes_rankproduct.tsv', sep='\t')
    return results

if __name__ == '__main__':
    # get_training_table()
    get_training_table(best_parameter_path='../ref_data/best_parameters_grid.csv', out_name='../genelist/Cancer_train_optimized.xlsx')
