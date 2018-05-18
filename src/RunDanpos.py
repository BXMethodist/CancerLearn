#!/usr/bin/env python
import os
from utils import *

if __name__ == '__main__':
    files = [x for x in os.listdir('/archive/tmhbxx3/CancerLearn/wigs/CD4/') if x.endswith('.wig')]

    for f in files:

        ## To generate wig
        # cmd = 'python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py dpeak ' + '/archive/tmhdxz9/work/h3k27me3/keji/wig/' + \
        #       f + ' -c 50000000 -f 0 --smooth_width 0 -o ' + f[:-4]

        ## To call peaks
        cmd = 'python /archive/tmhbxx3/tools/danpos_dev/danpos.py dpeak ' + '/archive/tmhbxx3/CancerLearn/wigs/CD4/' + \
              f + ' -c 50000000 -f 0 --smooth_width 0 -o ' + '/archive/tmhbxx3/CancerLearn/peaks/CD4/'+f[:-4]+' -q '
        cmd = cmd + ','.join([str(x) for x in range(1,101)])
        submit_pbs([cmd], f[:-4])

