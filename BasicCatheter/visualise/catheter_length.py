""" Module to make graphs from multiple catheter simulations. """

import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import seaborn as sns
import os
from matplotlib.font_manager import FontProperties

# Graphics settings
plt.rcParams["figure.figsize"] = (1.55,1.55)
plt.rcParams['xtick.major.pad']='4'
plt.rcParams['ytick.major.pad']='4'
font = FontProperties()
font.set_name('serif')
font.set_size(8)
colours = np.array([(68,120,33), (141,211,95), (255,179,128), (255,221,85), (179,179,179)])/255
palette = sns.color_palette(colours)
sns.set_palette(palette)
sns.set_style("ticks")
plt.rcParams['font.size'] = 8.0
plt.rcParams['font.family'] = 'serif'

# Location of results files
folder = '../ExploringCatheter/results/catheter_lengths'
files = os.listdir(folder)

# interested in time taken to infection (at outflow) against catheter length
infection_threshold = np.array([1,1e2,1e4]) # Clinical definition of infection occuring

# Arrays for results
N = len(files)
lengths = np.zeros(N)
times = np.zeros((N,len(infection_threshold)))

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
    r=float(info['outside growth rate'])
    D=float(info['surface diffusivity'])
    inverse_speed = 0.5/(np.sqrt(r*D)*3600*24)

    outflow = np.zeros(t_len)
    for j in range(0, t_len):
        outflow[j] = df.iloc[4*j+3][1]

    # Find times at which outflow crosses threshold
    indices = np.searchsorted(outflow,infection_threshold, side='right')
    times[i] = [dt*index/24 if index!=len(outflow) else np.nan for index in indices]

print(times)
print(inverse_speed)
#plt.plot(lengths, times)
plt.scatter(lengths,[time[0] for time in times], label='10$^3$ mm$^{{-3}}$',s=4,color=palette[0],zorder=2)
#plt.scatter(lengths,[time[1] for time in times], label='10$^5$ mm$^{{-3}}$',s=4,color=palette[1])
#plt.scatter(lengths,[time[2] for time in times], label='10$^7$ mm$^{{-3}}$',s=4,color=palette[3])
plt.plot(lengths, lengths*inverse_speed, color=palette[2], ls='--',zorder=1)
sns.despine()
plt.locator_params(tight=True, nbins=3)
plt.xticks(fontproperties=font)
plt.yticks(fontproperties=font)
plt.xlabel('Urethral length (mm)', fontproperties=font, labelpad=2, position=(0.475,0))
plt.ylabel('Time to detection (days)', fontproperties=font, labelpad=2, position=(0,0.435))
plt.tight_layout(rect=[-0.072,-0.068,1.075,1.07])
plt.savefig('catheter_lengths_timescale.pdf')
plt.show()