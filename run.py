#!/usr/bin/env python 

import SpotAnaliser
import sys
import os
import matplotlib.pyplot as plt
import prettyplotlib as ppl
import numpy as np

if len(sys.argv) < 3:
	print """USAGE: run.py PATH_TO_IMAGE_FILES TXT_IMAGE_FILE_LIST
	to run the example:
	python run.py example/images list.txt"""

path = os.path.join(os.getcwd(),sys.argv[1])
filelist = sys.argv[2]

# path = '/home/xqua/Documents/Harvard/Garner_Rotation/BMD203_Long_aquisition'
# lines = open('list2').readlines()

os.chdir(path)

lines = open(filelist).readlines()
bacterial_files = [i.strip() for i in lines if (i and 'tif' in i.lower())]

# filename = bacterial_files[25]
# filename = 'bMD369_NOiptg_1_R3D_nodrift-1.tif'

MreB_mu = []
MreB_std = []
lenght = []
speed = []

for filename in bacterial_files:
	A = SpotAnaliser.Spot_Analysis(filename)
	A.Main()
	if not np.isnan(A.meanMREB):
		MreB_mu.append(A.meanMREB)
		MreB_std.append(A.stdMREB)
		speed.append(A.speed)
		lenght.append(A.bact_lenght_nm)

fig = plt.figure(figsize=(15,15))
plt.errorbar(lenght, MreB_mu, yerr=MreB_std, fmt='o', alpha=0.7, color='g')
plt.xlabel('Bacteria Lenght in nm')
plt.ylabel('Number of detected MreB filaments')
plt.title('Correlation between the lenght of the bacteria and the number of MreB filaments')
fig.savefig('lenght_Filament_correlation.png')

