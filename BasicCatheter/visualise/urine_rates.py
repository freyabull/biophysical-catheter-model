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
folder = '../ExploringCatheter/results/urine_rates'
files = os.listdir(folder)

infection_threshold = np.array([1e3,1e5,1e7]) # Clinical definition of infection occuring

#  Threshold value for blockage
blockage_thresh = 1e6


# Arrays for results
N = len(files)
urine_rates = np.zeros(N)
times = np.zeros((N,len(infection_threshold)))
max_infect = np.zeros(N)
time_to_max = np.zeros(N)
time_to_block = np.zeros(N)

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

    bladder = np.zeros(t_len)
    inside = np.zeros((t_len,df.shape[1]))
    outflow = np.zeros(t_len)
    maxy = np.zeros(t_len)
    for j in range(0, t_len):
        bladder[j] = df.iloc[4*j+1][1]
        inside[j] = df.iloc[4*j+2]
        outflow[j] = df.iloc[4*j+3][1]
        maxy[j] = max(inside[j])

    # Find times at which outflow crosses threshold
    indices = np.searchsorted(outflow,infection_threshold, side='right')
    times[i] = [dt*index if index!=len(outflow) else np.nan for index in indices]

    # Find maximum bladder infection & time to infection
    max_infect[i] = np.max(bladder)
    time_to_max[i] = dt*np.searchsorted(bladder,0.99*max_infect[i], side='right')

    # Find times at which inside blocks
    index = np.searchsorted(maxy,blockage_thresh, side='right')
    time_to_block[i] = dt*index if index!=t_len else np.nan


fig, main_ax = plt.subplots()

main_ax.scatter(urine_rates,max_infect,zorder=2,color=palette[1])
main_ax.plot(urine_rates,max_infect,zorder=2)
sns.despine()
main_ax.set_yscale('log')
main_ax.axvline(50/3, color=palette[2], ls='--',zorder=1, label='Typical\npatient\nflow rate')
main_ax.legend(prop=font, bbox_to_anchor=(0.34, 0.21),handlelength=1,handletextpad=0.5,frameon=False)
plt.locator_params(axis='y',tight=True, numticks=4)
plt.locator_params(axis='x',tight=True, nbins=4)
plt.xticks(fontproperties=font)
plt.yticks(fontproperties=font)
main_ax.set_xlabel("Urine flow rate (mm$^3$s$^{-1}$)", fontproperties=font, labelpad=2)
main_ax.set_ylabel("Max bladder density (mm$^{-3}$)", fontproperties=font, labelpad=2)
sub_ax = fig.add_axes([0.72, 0.72, 0.25, 0.25])
plt.xticks(fontproperties=font)
plt.yticks(fontproperties=font)
sub_ax.plot(urine_rates, time_to_max,zorder=2)
sub_ax.axvline(50/3, color=palette[2], ls='--',zorder=1)
sub_ax.set_xlabel("Flow rate (mm$^3$s$^{-1}$)", fontproperties=font, labelpad=2,position=(0.355,0))
sub_ax.set_ylabel("Time to steady \n state (hr)", fontproperties=font, labelpad=2,position=(0,0.48))
plt.tight_layout(rect=[-0.045,-0.05,1.035,1.035])
plt.savefig('urine_rate_phases.pdf')

plt.figure()
plt.scatter(urine_rates,time_to_block)
sns.despine()
plt.locator_params(tight=True, nbins=4)
plt.xticks(fontproperties=font)
plt.yticks(fontproperties=font)
plt.xlabel("Urine flow rate (mm$^3$s$^{-1}$)", fontproperties=font, labelpad=2)
plt.ylabel('Time taken to internal blockage (hr)', fontproperties=font, labelpad=2)
plt.tight_layout(rect=[-0.04,-0.06,1.035,1.035])

plt.figure()
plt.plot(urine_rates, times) # Graph timescale to clinical infection at outlet against urine rate
plt.scatter(urine_rates,[time[0] for time in times], label='10$^3$ mm$^{{-3}}$')
plt.scatter(urine_rates,[time[1] for time in times], label='10$^5$ mm$^{{-3}}$')
plt.scatter(urine_rates,[time[2] for time in times], label='10$^7$ mm$^{{-3}}$')
plt.xticks(fontproperties=font)
plt.yticks(fontproperties=font)
plt.xlabel("Urine flow rate (mm s$^{-1}$", fontproperties=font, labelpad=2)
plt.ylabel("Time to detection (hr)", fontproperties=font, labelpad=2)
plt.legend(prop=font, title='Bacterial density')
legend_title_props = plt.gca().get_legend().get_title()
legend_title_props.set_name('serif')
legend_title_props.set_size(8)
plt.tight_layout(rect=[-0.04,-0.04,1.035,1.035])

plt.show()