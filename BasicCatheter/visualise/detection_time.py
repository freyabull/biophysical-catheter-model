""" Module to make graphs from multiple catheter simulations. """

import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import seaborn as sns
import os
from matplotlib.font_manager import FontProperties

# Graphics settings
plt.rcParams["figure.figsize"] = (3.2,1.6)
plt.rcParams['xtick.major.pad']='4'
plt.rcParams['ytick.major.pad']='4'
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = '8'
font = FontProperties()
font.set_name('serif')
font.set_size(8)
colours = np.array([(68,120,33), (141,211,95), (255,179,128), (255,221,85), (179,179,179)])/255
palette = sns.color_palette(colours)
palette2 = sns.dark_palette(colours[1], n_colors=4)
sns.set_palette(palette)
sns.set_style("ticks")
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = '8'

# Catheter lengths data processing

# Location of results files
folder = '../ExploringCatheter/results/catheter_lengths'
files = os.listdir(folder)

# interested in time taken to infection (at outflow) against catheter length
infection_threshold = np.array([1,1e2,1e4]) # Clinical definition of infection occuring

# Arrays for results
N = len(files)-2
lengths = np.zeros(N)
length_times = np.zeros((N,len(infection_threshold)))

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
    L=float(info['catheter length'])
    fisher_time = 0.5/(np.sqrt(r*D)*3600*24)

    outflow = np.zeros(t_len)
    for j in range(0, t_len):
        outflow[j] = df.iloc[4*j+3][1]

    # Find times at which outflow crosses threshold
    indices = np.searchsorted(outflow,infection_threshold, side='right')
    length_times[i] = [dt*index/24 if index!=len(outflow) else np.nan for index in indices]

# Urine rates processing

# Location of results files
folder = '../ExploringCatheter/results/urine_rates'
full_files = os.listdir(folder)
files = full_files[0:12:2]+full_files[12:16] + full_files[16:22:2] + full_files[22:]
print(files)

infection_threshold = np.array([1e0,1e2,1e4]) # Clinical definition of infection occuring

# Arrays for results
N = len(files)
urine_rates = np.zeros(N)
urine_times = np.zeros((N,len(infection_threshold)))

# Open files and read in data
for i in range(N):
    filename = os.path.join(folder, files[i])
    print(filename)
    # Read data
    df = pd.read_csv(filename, index_col=0, header=None, skiprows=3)
    info = pd.read_csv(filename,  nrows=1)
    if df.empty: break

    urine_rates[i] = info['urine rate']
    t_len = int(info['simulation length']/info['print interval']) # number of time steps
    dt = float(info['print interval']/3600) # print interval (hrs)

    outflow = np.zeros(t_len)
    for j in range(0, t_len):
        outflow[j] = df.iloc[4*j+3][1]

    # Find times at which outflow crosses threshold
    indices = np.searchsorted(outflow,infection_threshold, side='right')
    urine_times[i] = [dt*index if index!=len(outflow) else np.nan for index in indices]

# Set up figure/axes
gs_kw = dict(width_ratios=[1,1], height_ratios=[1])
big_fig, ax = plt.subplots(1,2, gridspec_kw=gs_kw)


ax[0].scatter(urine_rates,[time[0]/24 for time in urine_times], label='10$^3$ mm$^{{-3}}$',s=4,color=palette2[1],zorder=2)
ax[0].scatter(urine_rates,[time[1]/24 for time in urine_times], label='10$^5$ mm$^{{-3}}$',s=4,color=palette[0],zorder=2)
ax[0].scatter(urine_rates,[time[2]/24 for time in urine_times], label='10$^{7}$ mm$^{{-3}}$',s=4,color=palette[1],zorder=2)
ax[0].axhline(fisher_time*lengths[0], color=palette[2], ls='--',zorder=1)
ax[1].plot(lengths, lengths*fisher_time, color=palette[2], ls='--',zorder=1)
ax[1].scatter(lengths,[time[0] for time in length_times], label='10$^3$ mm$^{{-3}}$',s=4,color=palette[0],zorder=2)

sns.despine()
ax[0].locator_params(tight=True, nbins=3)
ax[1].locator_params(tight=True, nbins=3)
ax[0].set_xlabel('Urine rate (mm$^3$s$^{-1}$)', labelpad=11, va='bottom')
ax[1].set_xlabel('Urethral length (mm)', labelpad=11, va='bottom')
ax[0].set_ylabel('Detection time (days)', labelpad=2.0, fontproperties=font)
ax[1].set_ylabel('Detection time (days)', labelpad=2.0, fontproperties=font)
ax[1].set_xticks([60,120])

big_fig.add_subplot(111,frameon=False)
plt.tick_params(labelcolor='none',which='both',top=False,bottom=False,left=False,right=False)
plt.text(0.01,0.93, '(a)')
plt.text(0.61,0.93, '(b)')
plt.tight_layout(rect=[-0.115,-0.198,1.056,1.083])

plt.savefig('detection_time.pdf')