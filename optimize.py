#!/usr/bin/env python

import SpotAnaliser
from multiprocessing import Pool
import numpy as np
import scipy.optimize


def Loss_Function(res, known):
    return np.sum((res - known) ** 2)


def Fit(params, path, known):
    p = Pool(len(path))
    opt = []
    for i in path:
        opt.append([i, params])
    res = p.map(run, opt)
    L = Loss_Function(res, known)
    return L


def run(parameters):
    path, params = parameters
    n = int(path.split('/')[-1].split('_')[5])
    print "Nb of filament is ", n
    S = SpotAnaliser.Spot_Analysis(insilicopath + path)
    S.timethreshold = params[0]
    S.speed_threshold = params[1]
    S.Time_Error = params[2]
    S.Main()
    return S.totalMREB

if __name__ == '__main__':
    # insilicopath = "/n/Eggan_Lab3/Users/lblondel/"
    # path = ["sim_MreB_brightness_4000_numfilaments_40_speed_0.02_flen_0.2.mat",
    #         "sim_MreB_brightness_4000_numfilaments_10_speed_0.02_flen_0.2.mat",
    #         "sim_MreB_brightness_4000_numfilaments_80_speed_0.02_flen_0.2.mat",
    #         "sim_MreB_brightness_4000_numfilaments_120_speed_0.02_flen_0.2.mat",
    #         "sim_MreB_brightness_4000_numfilaments_180_speed_0.02_flen_0.2.mat",
    #         "sim_MreB_brightness_4000_numfilaments_10_speed_0.02_flen_0.5.mat",
    #         "sim_MreB_brightness_4000_numfilaments_40_speed_0.04_flen_0.2.mat",
    #         "sim_MreB_brightness_10000_numfilaments_80_speed_0.02_flen_0.2.mat",
    #         "sim_MreB_brightness_6000_numfilaments_120_speed_0.02_flen_0.2.mat",
    #         "sim_MreB_brightness_1000_numfilaments_40_speed_0.02_flen_0.2.mat"]
    insilicopath = "./"
    path = ["sim_MreB_brightness_4000_numfilaments_10_speed_0.01_flen_0.3.mat",
            "sim_MreB_brightness_4000_numfilaments_10_speed_0.01_flen_0.3.mat"]
    known = []
    for p in path:
        known.append(int(p.split('_')[5]))
    paths = [insilicopath + p for p in path]

    p0 = [60, 40, 50]
    bounds = [[0, 100], [0, 100], [0, 100]]
    scipy.optimize.minimize(Fit, p0, args=(paths, known), method="SLSQP", bounds=bounds, options={'maxiter': 500})
