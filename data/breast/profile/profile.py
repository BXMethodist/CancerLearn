import os, sys, pandas as pd

sys.path.append('/archive/tmhbxx3/tools/tools/')
from utils import *
from selector_utils import *
from scipy import stats

def profile(wigs, genefile_paths, wigfile_aliases, genefile_aliases, name, genomic_sites='TSS', flank_up=3000, flank_dn=3000, plot_colors=['blue','red','black'], heatmap=1):
    cmd = 'python /archive/tmhbxx3/tools/danpos_dev/danpos.py profile ' + ','.join(wigs)+' --genefile_paths '+','.join(genefile_paths)+' --wigfile_aliases '+','.join(wigfile_aliases) +' --genefile_aliases '+','.join(genefile_aliases) +' --name '+ name + ' --genomic_sites '+genomic_sites + ' --heatmap '+str(heatmap) +' --flank_up '+str(flank_up) + ' --flank_dn '+str(flank_dn) + ' --plot_colors '+','.join(plot_colors)
    return [cmd]

if __name__ == "__main__":
    """
    Get the profile of putative OGs and putative TSGs within normal breast cancer cells and tumor cancer cells
    """
    if False:
        normals = ['76NF2V', 'MCF10A']
        tumors = ['HCC1937', 'MB231',  'MB436', 'SKBR3', 'ZR751', 'AU565', 'HCC1954', 'MB361', 'MB468',  'MCF7','UACC812']
        # markers = ['H3K27ac', 'H3K4me3',  'H3K27me3',  'H3K4me1', 'H3K79me2',  'H3K9me3']
        markers = ['H3K27ac', 'H3K4me3',  'H3K27me3',  'H3K4me1']

        for marker in markers:
            for normal in normals:
                normal_wig = ['../'+marker+'/'+normal+'/'+x for x in os.listdir('../'+marker+'/'+normal+'/') if x.endswith('.wig')]
                for tumor in tumors:
                    tumor_wig = ['../'+marker+'/'+tumor+'/'+y for y in os.listdir('../'+marker+'/'+tumor+'/') if y.endswith('.wig')]
                    wigs = normal_wig + tumor_wig
                    cmds = profile(wigs, genefile_paths=[tumor+'_OGs.txt', tumor+'_TSGs.txt', tumor+'_controls.txt'], wigfile_aliases=[normal, tumor], genefile_aliases=['putative_OGs', 'putative_TSGs', 'controls'], name=os.getcwd()+'/'+normal+'_'+tumor+'_'+marker)
                    submit_pbs(cmds, normal+'_'+tumor+'_'+marker, mem="16000mb", queue='mediummem', ppn=1, walltime="32:00:00")

    """
    Calculate the profiles p values
    """
    if False:
        normals = ['76NF2V', 'MCF10A']
        tumors = ['HCC1937', 'MB231', 'MB436', 'SKBR3', 'ZR751', 'AU565', 'HCC1954', 'MB361', 'MB468', 'MCF7',
                  'UACC812']
        markers = ['H3K27ac', 'H3K4me3', 'H3K27me3', 'H3K4me1']

        results = []
        for normal in normals:
            for tumor in tumors:
                for marker in markers:
                    name = normal+'_'+tumor+'_'+marker+'_TSS.xls'
                    cur_df = pd.read_csv(name, sep='\t', index_col=0)
                    og_p = stats.mannwhitneyu(cur_df[tumor+'.putative_OGs.tss'], cur_df[normal+'.putative_OGs.tss'], alternative="greater")[1]
                    tsg_p = stats.mannwhitneyu(cur_df[tumor+'.putative_TSGs.tss'], cur_df[normal+'.putative_TSGs.tss'], alternative="less")[1]
                    control_p = stats.mannwhitneyu(cur_df[tumor+'.controls.tss'], cur_df[normal+'.controls.tss'])[1]
                    results.append([normal, tumor, marker, og_p, tsg_p, control_p])
        final_df = pd.DataFrame(results)
        final_df.columns = ['normal', 'tumor', 'marker', 'og_p', 'tsg_p', 'control_p']
        final_df.to_csv('profile_pvalues.xls', sep='\t', index=None)









