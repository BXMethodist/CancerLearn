import pandas as pd, os, numpy as np

if __name__ == "__main__":
    folders = [x for x in os.listdir('.') if os.path.isdir(x)]
    for folder in folders:
        if folder.find('Input') != -1:
            continue
        celltype, marker = folder.split('.')
        if not os.path.isdir('../wigs/'+marker):
            os.system('mkdir '+'../wigs/'+marker)
        if not os.path.isdir('../wigs/' + marker+'/'+celltype):
            os.system('mkdir '+'../wigs/'+marker+'/'+celltype)
        os.system('mv '+ folder+'/pooled/*.wig ../wigs/'+marker+'/'+celltype)
