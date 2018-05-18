import pandas as pd, os
from collections import defaultdict

if __name__ == "__main__":
    """
    combine the medata information and move the bowtie files into corresponding folders
    """
    if True:
        data_info = pd.read_csv('data_info.csv', index_col=0)
        meta_df = pd.read_csv('metadata.txt', sep='\t')
        meta_df = meta_df.set_index(['GSM_ID'])

        samples = defaultdict(dict)
        inputs = {}

        for GSM in meta_df.index:
            cur_data = data_info[data_info['SampleName']==GSM]
            if cur_data.shape[0]>1:
                cmd = 'cat '+','.join([x+'.bowtie' for x in cur_data.index])+' > '+ '_'.join(list(cur_data.index))+'.bowtie'
                move = '_'.join(list(cur_data.index))+'.bowtie'
            else:
                move = cur_data.index[0]+'.bowtie'
            # print meta_df.ix[GSM, 'information'].split(',')
            celltype, marker, _ = meta_df.ix[GSM, 'information'].split('.')
            if not os.path.isdir(celltype+'.'+marker):
                os.system('mkdir '+celltype+'.'+marker)
            os.system('cp '+move +' '+celltype+'.'+marker)
            if marker == 'Input':
                inputs[celltype] = celltype+'.'+marker
            else:
                samples[celltype][marker] = celltype+'.'+marker


