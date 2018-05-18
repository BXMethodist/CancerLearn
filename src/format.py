import os, pandas as pd

def change_name(path, before, after, **kwargs):
    if 'surfix' in kwargs.keys():
        surfix = kwargs['surfix']
    if surfix is not None:
        files = [x for x in os.listdir(path) if x.endswith(surfix)]
    else:
        return
    for f in files:
        os.system('mv '+path+'/'+f+' '+path+'/'+f.replace(before, after))


if __name__ == '__main__':
    change_name('/archive/tmhbxx3/CancerLearn/wigs/CD4', 'archive_tmhdxz9_work_h3k27me3_keji_wig_', 'keji_', surfix='.wig')
