#!/usr/bin/env python

import SpotAnaliser
import sys
import os
import numpy as np
import scipy.optimize

if len(sys.argv) < 3:
    print """USAGE: train.py PATH_TO_IMAGE_FILES TXT_IMAGE_FILE_LIST
    to run the example:
    python train.py example/images list_train.txt
    the file list_train.txt needs to have for each line: image_name,number of filament"""
    sys.exit()

path = os.path.join(os.getcwd(), sys.argv[1])
filelist = sys.argv[2]


os.chdir(path)

lines = open(filelist).readlines()
bacterial_files = [i.strip().split(',')[0]
                   for i in lines if (i and 'tif' in i.lower())]
known_filament_number = [float(i.strip().split(',')[1])
                         for i in lines if (i and 'tif' in i.lower())]


def loss_function(MreB, MreB_True):
    E = 0
    for i in range(len(MreB)):
        E += ((MreB[i] - MreB_True[i]) ** 2) / len(MreB)
    return E

timeinterval = 1.0
deltax = 0.1

timethreshold_totest = np.arange(4, 8, 0.5) * timeinterval * (1 / deltax)
timererror_totest = np.arange(40, 70, 10) / timeinterval
y_threshold_totest = np.arange(0.01, 0.1, 0.01)
dy_threshold_totest = np.arange(0, 0.0001, 0.00001)

Error = np.zeros((len(timethreshold_totest), len(timererror_totest), len(y_threshold_totest), len(dy_threshold_totest)))

for i in range(len(timethreshold_totest)):
    for j in range(len(timererror_totest)):
        for k in range(len(y_threshold_totest)):
            for l in range(len(dy_threshold_totest)):

                MreB_mu = []
                MreB_std = []
                lenght = []
                speed = []

                for filename in bacterial_files:
                    A = SpotAnaliser.Spot_Analysis(filename)
                    A.timethreshold = timethreshold_totest[i]
                    A.Time_Error = timererror_totest[j]
                    A.y_threshold = y_threshold_totest[k]
                    A.dy_threshold = dy_threshold_totest[l]

                    A.Main()
                    MreB_mu.append(A.meanMREB)
                    # MreB_std.append(A.stdMREB)
                    # speed.append(A.speed)
                    # lenght.append(A.bact_lenght_nm)
                E = loss_function(MreB_mu, known_filament_number)
                Error[i][j][k][l] = E

# print Error
# print np.where(Error == Error.min())
tmp = np.where(Error == Error.min())
i, j, k, l = tmp[0][0], tmp[1][0], tmp[2][0], tmp[3][0]
print "timethreshold = ", (timethreshold_totest[i] / timeinterval) / (1 / deltax)
print "Time_Error = ", timererror_totest[j] * timeinterval
print "y_threshold = ", y_threshold_totest[k]
print "dy_threshold = ", dy_threshold_totest[l]
