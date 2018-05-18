import os, sys


def submit_pbs(cmds, outputname):
    """
    :param cmd: command line
    :param outputname: output folder name
    :return:
    """
    pbs = open(outputname + ".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N " + outputname + '\n')
    pbs.write("#PBS -q mediummem\n")
    pbs.write("#PBS -m e\n")
    pbs.write("#PBS -M bxia@houstonmethodist.org\n")
    pbs.write("#PBS -l walltime=96:00:00\n")
    pbs.write("#PBS -l nodes=1:ppn=1\n")
    pbs.write("#PBS -l pmem=8000mb\n")
    pbs.write("cd " + os.getcwd() + "\n")
    pbs.write("module load python/2.7.11\n")
    for cmd in cmds:
        pbs.write(cmd + "\n")
    pbs.close()
    os.system('qsub ' + outputname + ".pbs")


if __name__ == "__main__":
    path = '/archive/tmhbxx3/CancerLearn/wigs/CD4/'
    wigs = [x for x in os.listdir(path) if x.endswith('.wig')]

    for i in range(0, len(wigs)):
        out_folder = '/archive/tmhbxx3/CancerLearn/peaks/CD4_sk/'+wigs[i][:-4].replace('keji_', '').replace('CD4-', '')
        os.system('mkdir '+out_folder)
        cmd = 'python /archive/tmhbxx3/CIG/test/skewness_kurtosis_run.py '+str(i)+' '+path+' '+ out_folder
        submit_pbs([cmd], wigs[i][:-4])
        #break