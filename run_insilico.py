#!/usr/bin/env python

import SpotAnaliser
import os
from multiprocessing import Pool
import time
import sys
import numpy as np

# l = os.listdir('InSilico')
insilicopath = "/n/Eggan_Lab3/Users/lblondel/"
respath = '/n/home05/lblondel/'
l = os.listdir(insilicopath)

l = [i for i in l if '.mat' in i and '.png' not in i and '.mp4' not in i]

done = []
try:
    f = open(respath + 'res_test', 'r')
    lines = f.readlines()
    for line in lines:
        s = line.strip().split('\t')
        done.append(s[0])
    f.close()
except:
    f = open(respath + 'res_test', 'w')
    legend = ['Path',
              'Real nb Mreb',
              'Total MREB',
              'Mean Sigma',
              'Std Sigma',
              'Mean Amplitude',
              'Std Ammplitude',
              'Mean Spread',
              'Std Spread',
              'Speed']
    f.write('\t'.join(legend) + '\n')

path = [i for i in l if i not in done]


def run(path):
    n = int(path.split('_')[5])
    print "Nb of filament is ", n
    if n < 150:
        S = SpotAnaliser.Spot_Analysis(insilicopath + path)
        S.Main()
        print "WRITTING !"
        if not os.path.isfile('lock'):
            print "no Lock, Writing now !"
            lock = open(respath + 'lock', 'w')
            lock.write('1')
            lock.close()
        else:
            while os.path.isfile(respath + 'lock'):
                print "lock on, waiting ..."
                time.sleep(1)

        f = open(respath + 'res_test', 'a')
        result = [path,
                  str(S.meta['num_particle'][0][0]),
                  str(S.totalMREB),
                  str(np.mean(S.All_Sigmas)),
                  str(np.std(S.All_Sigmas)),
                  str(np.mean(S.All_Intensity)),
                  str(np.std(S.All_Intensity)),
                  str(np.mean(S.All_Spread)),
                  str(np.std(S.All_Spread)),
                  str(S.speed)]
        txt = '\t'.join(result) + '\n'
        # txt = path + '\t' + str(S.meta['num_particle'][0][0]) + '\t' + str(S.meanMREB) + '\n'
        f.write(txt)
        f.close()
        print "result written removing lock ..."
        os.remove(respath + 'lock')
        return txt
    return None

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "USAGE: ./run_insilico.py NB_CORE"
        sys.exit()
    else:
        p = Pool(int(sys.argv[1]))
        p.map(run, path)
