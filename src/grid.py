from utils import submit_pbs
import pandas as pd

if __name__ == "__main__":
    df = pd.read_csv('../ref_data/parameters.csv')
    markers = df['marker'].unique()
    features = df['feature'].unique()

    # tables = ['../genelist/final_control_genes.xlsx', '../genelist/final_tsg_genes.xlsx']
    tables = ['../genelist/final_tsg_genes.xlsx']

    for table in tables:
        for marker in markers:
            for feature in features:
                if feature.find('single') != -1:
                    continue
                outname = marker+'_'+feature+'_'+table.split('_')[1]+'_vs_og'
                if feature.find('genebody') == -1:
                    if feature.find('kurtosis') != -1 or feature.find('skewness') != -1:
                        # continue
                        cmd = 'python /archive/tmhbxx3/tools/danpos_dev/danpos.py grid --TTS_pos TSS --up_stream_grid start:-100000:1000:1000:10000:2:1000 --down_stream_grid end:0:101000:1000:10000:2:1000 --height_grid height:1:51:1:2:2:1 -g 11 -n gene -t cell_type ' + table + ' ../genelist/final_og_genes.xlsx /archive/tmhbxx3/CancerLearn/peaks/CD4_sk/' + marker + '/ ' + feature + ' ' + outname + ' ../grid/ /archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.xls'
                    else:
                        cmd = 'python /archive/tmhbxx3/tools/danpos_dev/danpos.py grid --TTS_pos TSS --up_stream_grid start:-100000:1000:1000:10000:2:1000 --down_stream_grid end:0:101000:1000:10000:2:1000 --height_grid height:1:51:1:2:2:1 -g 11 -n gene -t cell_type '+table+' ../genelist/final_og_genes.xlsx /archive/tmhbxx3/CancerLearn/peaks/CD4/'+marker+'/ '+feature+' '+outname+' ../grid/ /archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.xls'
                else:
                    feature = feature.replace('_genebody', '')
                    cmd = 'python /archive/tmhbxx3/tools/danpos_dev/danpos.py grid --TTS_pos TTS --up_stream_grid start:-50000:1000:1000:5000:2:1000 --down_stream_grid end:0:51000:1000:5000:2:1000 --height_grid height:1:51:1:2:2:1 -g 11 -n gene -t cell_type ' + table + ' ../genelist/final_control_genes.xlsx /archive/tmhbxx3/CancerLearn/peaks/CD4/' + marker + '/ ' + feature + ' ' + outname + ' ../grid/ /archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.xls'

                submit_pbs([cmd], outname, queue='default', time='96:00:00', ppn='8', mem="500mb")