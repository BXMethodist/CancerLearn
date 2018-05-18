import os
from utils import submit_pbs

if __name__ == '__main__':
    files = [x for x in os.listdir('/archive/tmhbxx3/CancerLearn/data/') if x.startswith('SRR')]

    cmd = []
    for f in files:
        cmd.append('/archive/tmhbxx3/tools/sratoolkit/bin/fastq-dump '+'/archive/tmhbxx3/CancerLearn/data/'+f)
    submit_pbs(cmd, 'cancer', queue='default', time='96:00:00', mem='4000mb')