""" Module to make graphs from multiple catheter simulations. """

import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import seaborn as sns
import os
from matplotlib.font_manager import FontProperties

# Graphics settings
#plt.rcParams["figure.figsize"] = (11/2.54,8/2.54)
plt.rcParams["figure.figsize"] = (3.1,3.1)
plt.rcParams['xtick.major.pad']='4'
plt.rcParams['ytick.major.pad']='4'
font = FontProperties()
#font.set_name('sans-serif')
font.set_name('serif')
font.set_size(8)
colours = np.array([(68,120,33), (141,211,95), (255,179,128), (255,221,85), (179,179,179)])/255
palette = sns.color_palette(colours)
sns.set_palette(palette)
sns.set_style("ticks")


# Location of results files
folder = '../ExploringCatheter/results/sump_volumes'
files = os.listdir(folder)

# Arrays for results
N = len(files)
sump_volumes = np.zeros(N)
times = np.zeros(N)

#  Threshold value for blockage
thresh = 1e6


# Open files and read in data
for i in range(N):
    filename = os.path.join(folder, files[i])
    print(filename)
    # Read data
    df = pd.read_csv(filename, index_col=0, header=None, skiprows=3)
    info = pd.read_csv(filename,  nrows=1)
    if df.empty: break

    sump_volumes[i] = info['sump volume']
    t_len = int(info['simulation length']/info['print interval']) # number of time steps
    dt = float(info['print interval']/3600) # print interval (hrs)

    inside = np.zeros((t_len,df.shape[1]))
    maxy = np.zeros(t_len)
    for j in range(0, t_len):
        inside[j] = df.iloc[4*j+2]
        maxy[j] = max(inside[j])

    index = np.searchsorted(maxy,thresh, side='right')
    times[i] = dt*index if index!=t_len else np.nan

plt.scatter(sump_volumes,times)
sns.despine()
plt.xscale('log')
plt.locator_params(axis='x',tight=True, numticks=4)
plt.locator_params(axis='y',tight=True, nbins=4)
plt.xticks(fontproperties=font)
plt.yticks(fontproperties=font)
plt.xlabel('Volume of residual urine (mm$^3$)', fontproperties=font, labelpad=2)
plt.ylabel('Time taken to internal blockage (hr)', fontproperties=font, labelpad=2)
plt.tight_layout(rect=[-0.04,-0.06,1.035,1.035])
plt.savefig('sump_volumes_blockages.pdf')
plt.show()

