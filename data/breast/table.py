"""
This module take a table with best parameter combinations, and generate the table with features for each genes.
"""
import os, sys, pandas as pd, numpy as np, random
from collections import defaultdict
from copy import deepcopy
from multiprocessing import Process, Queue
from itertools import permutations
from scipy.stats import norm

def random_control_genes(CIG_gene_df, exclude_list, all_genes, random_seed, number_genes=None, candidates=None):
    """
    random select genes for each cell type

    :param CIG_gene_table: table path
    :param exclude_list: a dict containing exclusion genes name for each cell type
    :param all_genes: all gene GTF
    :param random_seed: random seed number to make sure each return the same random gene set
    :return:
    """
    results = []
    random.seed(random_seed)
    for cell_type in CIG_gene_df['cell_type'].unique():
        if candidates is not None:
            cur_df = CIG_gene_df[CIG_gene_df['cell_type'] == cell_type]
            cur_number = cur_df.shape[0]
            if number_genes is not None:
                cur_number = int(cur_number*number_genes)
            else:
                cur_number = cur_number
            cur_negative_df = candidates[candidates.cell_type==cell_type].sample(cur_number, random_seed=random_seed)
            for j in range(cur_negative_df.shape[0]):
                results.append(cur_negative_df.iloc[j, :])
        else:
            cur_df = CIG_gene_df[CIG_gene_df['cell_type']==cell_type]
            cur_number = cur_df.shape[0]
            cur_candidates = set()
            for gene in all_genes:
                if isinstance(exclude_list, dict):
                    if gene not in exclude_list[cell_type] and gene not in CIG_gene_df['gene'].unique() and not gene.startswith('MIR'):
                        cur_candidates.add(gene)
                elif isinstance(exclude_list, set):
                    # print 'lala'
                    if gene not in exclude_list and gene not in CIG_gene_df['gene'].unique() and not gene.startswith('MIR'):
                        cur_candidates.add(gene)

            # print 'ISOC1' in cur_candidates, 'why CIG in?'

            if number_genes is not None:
                cur_negative = random.sample(cur_candidates, int(cur_number*number_genes))
            else:
                cur_negative = random.sample(cur_candidates, cur_number)
            for i in range(len(cur_negative)):
                try:
                    results.append([cur_negative[i]]+list(cur_df.iloc[i, 1:]))
                except:
                    results.append([cur_negative[i]] + list(cur_df.iloc[0, 1:]))


    # candidates = set()
    # # all_genes = set(pd.read_csv('hg19.GREATgene2UCSCknownGenes.table.xls', sep='\t', index_col=0).index)
    # for gene in all_genes:
    #     if gene not in exclude_list and gene not in CIG_gene_df['gene'] and not gene.startswith('MIR'):
    #         candidates.add(gene)
    # candidates = list(candidates)
    # random.seed(random_seed)
    #
    # if number_genes is None:
    #     number_genes = len(candidates)
    #
    # negative_control_genes = random.sample(candidates, number_genes)
    #
    # for i in range(len(negative_control_genes)):
    #     try:
    #         results.append([negative_control_genes[i]]+list(CIG_gene_df.iloc[i, 1:]))
    #     except:
    #         results.append([negative_control_genes[i]] + list(CIG_gene_df.iloc[0, 1:]))
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

def get_stats(gene_df, df_path, criteria, cell_type, bin=3000, df_function=df_to_index_danpos):
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
        final.append((gene_name, results[gene_name], cell_type))
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

    for cell_type in CIG_gene_df['cell_type'].unique():
        cur_CIG_gene_list = list(CIG_gene_df[CIG_gene_df['cell_type'] == cell_type]['gene'].values)
        cur_CIG_gene_range = CIG_gene_ranges[CIG_gene_ranges['gene'].isin(cur_CIG_gene_list)]

        cur_non_CIG_gene_list = list(non_CIG_gene_df[non_CIG_gene_df['cell_type'] == cell_type]['gene'].values)
        cur_non_CIG_gene_range = non_CIG_gene_ranges[non_CIG_gene_ranges['gene'].isin(cur_non_CIG_gene_list)]

        cur_df = all_dfs[cell_type][cutoff]

        cur_CIG_result = get_stats(cur_CIG_gene_range, cur_df, criteria, cell_type)
        cur_non_CIG_result = get_stats(cur_non_CIG_gene_range, cur_df, criteria, cell_type)

        CIG_results += cur_CIG_result
        non_CIG_results += cur_non_CIG_result
    CIG_results_df = pd.DataFrame(CIG_results)

    if TTS_pos == 'TTS':
        CIG_results_df.columns = ['gene', marker + '_' + criteria + '_' + 'genebody', 'cell_type']
    else:
        CIG_results_df.columns = ['gene', marker + '_' + criteria, 'cell_type']
    CIG_results_df['name'] = CIG_results_df['gene'] + '_'+ CIG_results_df.cell_type
    CIG_results_df=CIG_results_df.set_index(['name'])

    non_CIG_results_df = pd.DataFrame(non_CIG_results)
    if TTS_pos == 'TTS':
        non_CIG_results_df.columns = ['gene', marker + '_' + criteria + '_' + 'genebody', 'cell_type']
    else:
        non_CIG_results_df.columns = ['gene', marker +'_'+ criteria, 'cell_type']
    non_CIG_results_df['name'] = non_CIG_results_df['gene'] + '_' + non_CIG_results_df.cell_type
    non_CIG_results_df = non_CIG_results_df.set_index(['name'])

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

def CIG_exclude_list(CIG_gene_df, CIG_exclude_genes=None):
    if CIG_exclude_genes is None:
        return CIG_gene_df
    CIG_gene_df = CIG_gene_df[~CIG_gene_df['gene'].isin(CIG_exclude_genes)]
    return CIG_gene_df

def get_training_table(all_gene_GTF_path='/archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.xls',
                       exclude_list_file='/home/tmhbxx3/scratch/CIG/test/related_genes/',
                       CIG_gene_df_path='CIG_strong_bo.csv',
                       best_parameter_path='best_parameters_CIG_6th_bo_1xControl.csv', CIG_exclude_genes=None, candidates=None,
                       process=8, size=1, random_seed=1,
                       out_name=None):
    ## get the CIG genesets and nonCIG genesets
    all_gene_GTF = pd.read_csv(all_gene_GTF_path, sep='\t', dtype={'hg19.kgXref.geneSymbol': str})

    # all_gene_GTF = pd.read_csv('hg19.GREATgene2UCSCknownGenes.table.xls', sep='\t')
    all_gene_GTF['hg19.kgXref.geneSymbol'] = all_gene_GTF['hg19.kgXref.geneSymbol'].str.upper()
    # all_gene_GTF = all_gene_GTF[(all_gene_GTF['hg19.knownGene.txEnd'] - all_gene_GTF['hg19.knownGene.txStart'])>200]

    all_genes = set(all_gene_GTF['hg19.kgXref.geneSymbol'].unique())

    if exclude_list_file is not None:
        exclude_list = defaultdict(set)
        exclude_list_file = [exclude_list_file + x for x in os.listdir(exclude_list_file)]
        for related in exclude_list_file:
            related_file = open(related, 'r')
            related_type = related[related.rfind('/')+1:-4]
            related_list = [x.strip() for x in related_file]
            related_file.close()
            exclude_list[related_type] = set(related_list)
    else:
        f = open('merge_ten_cells_relative_genes.txt', 'r')
        exclude_list = set([g.strip().upper() for g in f.readlines()])
        f.close()
        # print exclude_list

    if CIG_gene_df_path.endswith('.csv'):
        CIG_gene_df = pd.read_csv(CIG_gene_df_path, dtype={'gene': str})
    elif CIG_gene_df_path.endswith('.xlsx'):
        CIG_gene_df = pd.read_excel(CIG_gene_df_path)
    # CIG_gene_df = pd.read_csv('top500.tuson.oncogene_keji.csv')
    # CIG_gene_df = pd.read_csv('top500.tuson.tumorSuppressor_keji.csv')
    CIG_gene_df['gene'] = CIG_gene_df['gene'].str.upper()
    CIG_gene_df = CIG_exclude_list(CIG_gene_df, CIG_exclude_genes)

    # print CIG_gene_df.shape, 'CIG genes'

    ## for grid optimization
    # non_CIG_gene_df = random_control_genes(CIG_gene_df, exclude_list, all_genes, random_seed=1, number_genes=size, candidates=candidates)
    # non_CIG_gene_df['gene'] = non_CIG_gene_df['gene'].str.upper()
    # non_CIG_gene_df.to_csv('non_CIG_control_v6_grid.csv', index=None, encoding='utf-8')

    ## for lg
    non_CIG_gene_df = random_control_genes(CIG_gene_df, exclude_list, all_genes, random_seed=random_seed, number_genes=size)
    non_CIG_gene_df['gene'] = non_CIG_gene_df['gene'].str.upper()
    non_CIG_gene_df['name'] = non_CIG_gene_df['gene'] + '_' + non_CIG_gene_df['cell_type']
    non_CIG_gene_df.to_csv('non_CIG_control_v311_size'+str(size)+'.csv', index=None, encoding='utf-8')

    # return
    # non_CIG_gene_df = pd.read_csv('non_CIG_control.csv', dtype={'gene': str})
    # non_CIG_gene_df['gene'] = non_CIG_gene_df['gene'].str.upper()

    # get feature parameter table
    best_feature_df = pd.read_csv(best_parameter_path)

    # get all tables:
    all_tables = {}
    all_sk_tables = {}

    for m in ['h3k4me3_qn', 'h3k27me3_qn', 'h3k4me1_qn', 'h3k27ac_qn']:
        if m not in all_tables:
            all_tables[m] = defaultdict(dict)
            dfs_path = '/home/tmhbxx3/scratch/CIG/' + m + '_peaks/pooled/'
            # if best_parameter_path.find('Great') == -1:
            #     dfs_path = '/home/tmhbxx3/scratch/CIG/' + m + '_peaks/pooled/'
            # else:
            #     dfs_path = '/home/tmhbxx3/scratch/CIG/' + m + 'default_peaks/pooled/'
            dfs = [x for x in os.listdir(dfs_path) if x.endswith('.xls')]

            for table_name in dfs:
                info = table_name.split('_')
                # print info
                cell_type = info[2]
                cutoff = float(info[-1][:-4])
                all_tables[m][cell_type][cutoff] = dfs_path+table_name

        if m not in all_sk_tables:
            all_sk_tables[m] = defaultdict(dict)
            dfs_path = '/home/tmhbxx3/scratch/CIG/test/' + m + '_sk_peaks/'
            dfs = [x for x in os.listdir(dfs_path) if x.endswith('.csv')]

            for table_name in dfs:
                info = table_name.split('_')
                cell_type = info[0]
                cutoff = float(info[-1][:-4])
                all_sk_tables[m][cell_type][cutoff] = dfs_path + table_name

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

    final_feature_df['new'] = final_feature_df.index
    two_columns = final_feature_df.new.apply(lambda x: pd.Series(x.split('_')))

    final_feature_df['cell_type'] = two_columns.iloc[:, 1]
    # print final_feature_df
    final_feature_df['gene_id'] = two_columns.iloc[:, 0]

    del final_feature_df['new']

    print final_feature_df.columns

    final_columns = [x for x in final_feature_df.columns if x.find('gene')== -1 or x == 'gene_id' or x.find('genebody')!=-1]

    final_feature_df = final_feature_df[final_columns]

    celltype_genes = final_feature_df['gene_id']+ final_feature_df['cell_type']

    cig_celltype_genes = CIG_gene_df['gene'] + CIG_gene_df['cell_type']
    for gene in cig_celltype_genes.unique():
        if gene not in celltype_genes.unique():
            print 'celltype wrong',  gene

    final_feature_df = final_feature_df.set_index(['gene_id'])
    if out_name is None:
        final_feature_df.to_csv('CIG_train_6th.csv')
    else:
        final_feature_df.to_csv(out_name)


def training_table_process(queue, best_feature_df, CIG_gene_df, non_CIG_gene_df, all_gene_GTF, all_tables, all_sk_tables):
    cur_final_feature_df = pd.DataFrame(index=list(CIG_gene_df.name.values) + list(non_CIG_gene_df.name.values))
    # print best_feature_df.columns
    for i in range(best_feature_df.shape[0]):
        marker = best_feature_df.iloc[i, 0] + '_qn'
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


def get_real_table(target_cell_type, all_gene_GTF_path='/archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.xls',
                   best_parameter_path='parameters_default.csv',
                   peak_path='/archive2/tmhbxx3/CancerLearn/data/normal/wigs_qn/peaks',
                   process=8,
                   out_name=None):
    ## get the CIG genesets and nonCIG genesets
    all_gene_GTF = pd.read_csv(all_gene_GTF_path, sep='\t', dtype={'hg19.kgXref.geneSymbol': str})

    # all_gene_GTF = pd.read_csv('hg19.GREATgene2UCSCknownGenes.table.xls', sep='\t')
    all_gene_GTF['hg19.kgXref.geneSymbol'] = all_gene_GTF['hg19.kgXref.geneSymbol'].str.upper()

    all_genes = set(all_gene_GTF['hg19.kgXref.geneSymbol'].unique())

    # get feature parameter table
    best_feature_df = pd.read_csv(best_parameter_path)

    # get all tables:
    all_tables = {}
    all_sk_tables = {}

    m_types = os.listdir(peak_path+'/'+target_cell_type)
    for m in m_types:
        if m not in all_tables:
            all_tables[m] = defaultdict(dict)
            dfs_path = peak_path+'/'+target_cell_type+'/'+m+'/'
            dfs = [x for x in os.listdir(dfs_path) if x.endswith('.xls')]

            for table_name in dfs:
                info = table_name.split('_')
                cell_type = target_cell_type
                cutoff = float(info[-1][:-4])
                all_tables[m][cell_type][cutoff] = dfs_path+table_name

        if m not in all_sk_tables:
            all_sk_tables[m] = defaultdict(dict)
            dfs_path = peak_path+'_sk/'+target_cell_type+'/'+m+'/'
            dfs = [x for x in os.listdir(dfs_path) if x.endswith('.csv')]

            for table_name in dfs:
                info = table_name.split('_')
                cell_type = target_cell_type
                cutoff = float(info[-1][:-4])
                all_sk_tables[m][cell_type][cutoff] = dfs_path + table_name

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

    queue = Queue()
    processes = []

    for i in range(process):
        cur_chunk = chunks[i]
        p = Process(target=real_table_process,
                    args=(queue, target_cell_type, cur_chunk, all_gene_GTF, all_tables, all_sk_tables, m_types))
        processes.append(p)
        p.start()

    final_feature_tables = []

    for i in range(process):
        cur_final_feature_df = queue.get()
        final_feature_tables.append(cur_final_feature_df)
    final_feature_df = pd.concat(final_feature_tables, axis=1)

    for p in processes:
        p.join()

    index_cols = []
    index_col = False
    for ci in range(len(final_feature_df.columns)):
        if final_feature_df.columns[ci].startswith('gene'):
            if not index_col:
                index_cols.append(ci)
                index_col = True
        else:
            index_cols.append(ci)
    final_feature_df = final_feature_df.iloc[:, index_cols]
    columns = list(final_feature_df.columns)
    columns[0] = 'gene_id'
    final_feature_df.columns = columns

    if out_name is None:
        final_feature_df.to_csv(target_cell_type+'_' + best_parameter_path[best_parameter_path.rfind('/')+1:-4]+'_'+'real_table.csv', index=None)
    else:
        final_feature_df.to_csv(out_name, index=None)

def real_table_process(queue, target_cell_type, best_feature_df, all_gene_GTF, all_tables, all_sk_tables, m_types):
    cur_final_feature_df = None
    # print best_feature_df.columns
    for i in range(best_feature_df.shape[0]):
        marker = best_feature_df.iloc[i, 0]
        if marker not in m_types:
            continue
        feature = best_feature_df.iloc[i, 1]

        criteria = feature.replace('_genebody', '')

        start = best_feature_df.iloc[i, 2]
        end = best_feature_df.iloc[i, 3]
        height = best_feature_df.iloc[i, 4]
        # print marker, feature, start, end, height
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

        all_gene_results_df = CIG_selecter_all(all_gene_GTF, start, end, cur_stat_dfs, height,
                                               criteria, target_cell_type, TSS, TTS)

        all_gene_results_df.columns = ['gene'+str(i), marker + '_' + feature]
        if cur_final_feature_df is None:
            cur_final_feature_df = all_gene_results_df.copy()
        else:
            cur_final_feature_df[marker + '_' + feature] = all_gene_results_df[all_gene_results_df.columns[1]]

    # cur_final_feature_df.to_csv(str(os.getpid())+'_'+target_cell_type+'_'+'real_table.csv')
    queue.put(cur_final_feature_df)
    return

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
    index = int(sys.argv[1])
    cell_types = os.listdir('/archive2/tmhbxx3/CancerLearn/data/breast/peaks/')
    cell_type = cell_types[index]

    get_real_table(cell_type,
                   peak_path='/archive2/tmhbxx3/CancerLearn/data/breast/peaks')



