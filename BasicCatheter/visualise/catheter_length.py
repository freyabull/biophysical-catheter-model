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
folder = '../ExploringCatheter/results/catheter_lengths'
files = os.listdir(folder)

# Arrays for results
N = len(files)
lengths = np.zeros(N)
times = np.zeros(N)

# interested in time taken to infection (at outflow) against catheter length
thresh = 1e5

for i in range(N):
    filename = os.path.join(folder, files[i])
    print(filename)
    # Read data
    df = pd.read_csv(filename, index_col=0, header=None, skiprows=3)
    info = pd.read_csv(filename,  nrows=1)
    if df.empty: break

    lengths[i] = info['catheter length']
    t_len = int(info['simulation length']/info['print interval']) # number of time steps
    dt = float(info['print interval']/3600) # print interval (hrs)

    outflow = np.zeros(t_len)
    for j in range(0, t_len):
        outflow[j] = df.iloc[4*j+3][1]

    index = np.searchsorted(outflow,thresh, side='right')
    times[i] = dt*index if index!=len(outflow) else np.nan


plt.scatter(lengths, times)
sns.despine()
plt.locator_params(tight=True, nbins=4)
plt.xticks(fontproperties=font)
plt.yticks(fontproperties=font)
plt.xlabel('Urethral length (mm)', fontproperties=font, labelpad=2)
plt.ylabel('Time till infection observed (hr)', fontproperties=font, labelpad=2)
plt.tight_layout(rect=[-0.04,-0.04,1.04,1.04])
plt.savefig('catheter_lengths_timescale.pdf')
plt.show()