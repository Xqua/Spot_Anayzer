#!/usr/bin/env python

import SpotAnaliser
import os
from multiprocessing import Pool
import time
import sys

l = os.listdir('InSilico')
l = [i for i in l if '.mat' in i and '.png' not in i and '.mp4' not in i]

done = []
try:
    f = open('res_test', 'r')
    lines = f.readlines()
    for line in lines:
        s = line.strip().split('\t')
        done.append(s[0])
    f.close()
except:
    pass

path = [i for i in l if i not in done]


def run(path):
    n = int(path.split('_')[-1].split('.')[0])
    if n < 150:
        S = SpotAnaliser.Spot_Analysis('InSilico/' + path)
        S.Main()
        print "WRITTING !"
        if not os.path.isfile('lock'):
            lock = open('lock', 'w')
            lock.write('1')
            lock.close()
        else:
            while os.path.isfile('lock'):
                time.sleep(1)

        f = open('res_test', 'a')
        txt = path + '\t' + str(S.meta['num_particle'][0][0]) + '\t' + str(S.meanMREB) + '\n'
        f.write(txt)
        f.close()

        os.remove('lock')
        return txt
    return None

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "USAGE: ./run_insilico.py NB_CORE"
        sys.exit()
    else:
        p = Pool(sys.argv[1])
        p.map(run, path)
