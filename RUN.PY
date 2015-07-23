#!/usr/bin/env python

import SpotAnaliser
import os
from multiprocessing import Pool
import time
import sys
import numpy as np

# l = os.listdir('InSilico')
# THIS IS THE DIRECTORY CONTAINING THE DIFFERENT MOVIES TO BE ANALYZED !
insilicopath = "/n/Eggan_Lab3/Users/lblondel/"

# THIS IS THE DIRECTORY WHERE THE RESULT WILL BE WRITTEN !
respath = '/n/home05/lblondel/'

# THIS IS TO READ THE DIRECTORY AND GET THE LIST OF THE FILES IN IT THEN FILTER OUT OF THERE IS JUNK
l = os.listdir(insilicopath)  # LIST DIR
l = [i for i in l if '.mat' in i and '.png' not in i and '.mp4' not in i]  # FILTER OUT JUNK

# THIS FUNCTION CHECKS FOR ALREADY DONE FILES SO THAT IT DOES NOT REDO ALREADY ANALYSED FILES.
# THIS IS FOR THE IN CASE OF IT CRASHES, RESTART FROM THE LAST ONE ANALYZED !
# OR IF THE RESULT FILE DOES NOT EXIST IT CREATES IT !
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
        S = SpotAnaliser.Spot_Analysis(insilicopath + path)
        S.Main()
        f = open(respath + 'res_test', 'a')
        result = [path,
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
        return txt
    return None

if __name__ == '__main__':
    for p in path:
        r = run(p)
        print r
