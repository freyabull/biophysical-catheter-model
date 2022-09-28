""" Module to make graphs from multiple catheter simulations. """

import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import seaborn as sns
import os
from matplotlib.font_manager import FontProperties

# Graphics settings
plt.rcParams["figure.figsize"] = (3.1,3.1)
plt.rcParams['xtick.major.pad']='4'
plt.rcParams['ytick.major.pad']='4'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = '8'
font = FontProperties()
font.set_name('serif')
font.set_size(8)
colours = np.array([(68,120,33), (141,211,95), (255,179,128), (255,221,85), (179,179,179)])/255
palette = sns.color_palette(colours)
sns.set_palette(palette)
sns.set_style("ticks")

# Location of results files
folder = '../ExploringCatheter/results/urine_rates'
files = os.listdir(folder)

infection_threshold = np.array([1e0,1e2,1e4]) # Clinical definition of infection occuring

#  Threshold value for blockage
blockage_thresh = 1e7


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
    r=float(info['outside growth rate'])
    D=float(info['surface diffusivity'])
    L=float(info['catheter length'])
    fisher_time = 0.5*L/(np.sqrt(r*D)*3600*24)

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
sub_ax = fig.add_axes([0.739, 0.745, 0.25, 0.25])
sub_ax2 = fig.add_axes([0.739, 0.47, 0.25, 0.25], sharex=sub_ax)
sub_ax.tick_params(labelbottom=False)
plt.xticks(fontproperties=font)
plt.yticks(fontproperties=font)
sub_ax.plot(urine_rates, time_to_max/24,zorder=2)
sub_ax.axvline(50/3, color=palette[2], ls='--',zorder=1)
sub_ax.set_yticks([170,190])
sub_ax.set_ylim((167,191))
sub_ax.set_ylabel("Time (days)", fontproperties=font, labelpad=2.5,position=(0,0.53))
sub_ax2.plot(urine_rates,max_infect)
sub_ax2.axvline(50/3, color=palette[2], ls='--',zorder=1)
sub_ax2.set_yticks([0,1e6])
sub_ax2.set_yticklabels(['0','$10^6$'])
sub_ax2.set_ylim((-5e4,105e4))
sub_ax2.set_ylabel("Density (mm$^{-3}$)", fontproperties=font, labelpad=0, position=(0,0.46))
sub_ax2.set_xlabel("Urine (mm$^3$s$^{-1}$)", fontproperties=font, labelpad=2,position=(0.48,0))
plt.tight_layout(rect=[-0.043,-0.044,1.033,1.036])
plt.savefig('urine_rate_phases.pdf')




plt.rcParams["figure.figsize"] = (1.55,1.55)

plt.figure()
plt.scatter(urine_rates,[time[0]/24 for time in times], label='10$^3$ mm$^{{-3}}$',s=4,color=palette[0],zorder=2)
plt.scatter(urine_rates,[time[1]/24 for time in times], label='10$^5$ mm$^{{-3}}$',s=4,color=palette[1],zorder=2)
plt.scatter(urine_rates,[time[2]/24 for time in times], label='10$^{7}$ mm$^{{-3}}$',s=4,color=palette[3],zorder=2)
plt.axhline(fisher_time, color=palette[2], ls='--',zorder=1)
sns.despine()
plt.locator_params(tight=True, nbins=3)
plt.xticks(fontproperties=font)
plt.yticks(fontproperties=font)
plt.xlabel("Urine rate (mm$^3$s$^{-1}$)", fontproperties=font, labelpad=1.5, position=(0.49,0))
plt.ylabel("Time to detection (days)", fontproperties=font, labelpad=2,position=(0,0.43))
plt.tight_layout(rect=[-0.071,-0.092,1.076,1.069]) #lbrt
plt.savefig('urine_rate_detection.pdf')


plt.figure()
plt.scatter(urine_rates,time_to_block,s=4,color=palette[0])
sns.despine()
plt.ylim(4033,4073)
plt.locator_params(tight=True, nbins=3)
plt.xticks(ticks=[20,40],fontproperties=font)
plt.yticks(ticks=[4040,4060],fontproperties=font)
plt.xlabel("Urine rate (mm$^3$s$^{-1}$)", fontproperties=font, labelpad=1.5, position=(0.47,0))
plt.ylabel('Time to thick biofilm (hr)', fontproperties=font, labelpad=2, position=(0,0.42))
plt.tight_layout(rect=[-0.072,-0.092,1.076,1.069])
plt.savefig('urine_rate_blockage.pdf')

#plt.rcParams["figure.figsize"] = (3.1,3.1)
#fig, main_ax = plt.subplots()

#main_ax.scatter(urine_rates,max_infect,zorder=2,color=palette[1])
#main_ax.plot(urine_rates,max_infect,zorder=2)
#sns.despine()
#main_ax.set_yscale('log')
#main_ax.axvline(50/3, color=palette[2], ls='--',zorder=1, label='Typical\npatient')
#main_ax.legend(prop=font, bbox_to_anchor=(0.34, 0.21),handlelength=1,handletextpad=0.5,frameon=False)
#plt.locator_params(axis='y',tight=True, numticks=4)
#plt.locator_params(axis='x',tight=True, nbins=4)
#plt.xticks(fontproperties=font)
#plt.yticks(fontproperties=font)
#main_ax.set_xlabel("Urine production rate (mm$^3$s$^{-1}$)", fontproperties=font, labelpad=2)
#main_ax.set_ylabel("Bacterial density in bladder (mm$^{-3}$)", fontproperties=font, labelpad=2)
#plt.xticks(fontproperties=font)
#plt.yticks(fontproperties=font)
#plt.tight_layout(rect=[-0.043,-0.044,1.033,1.036])
#plt.savefig('urine_rate_3mt.svg')


#plt.rcParams['font.size'] = '12'
#font2 = FontProperties()
#font2.set_name('serif')
#font2.set_size(12)

#fig, main_ax = plt.subplots()
#main_ax.scatter(urine_rates,max_infect,zorder=2,color=palette[1])
#main_ax.plot(urine_rates,max_infect,zorder=2)
#sns.despine()
#main_ax.set_yscale('log')
#main_ax.axvline(50/3, color=palette[2], ls='--',zorder=1, label='Typical\npatient')
#annie = main_ax.annotate('Typical patient', xy=(17,4e5), xytext=(22,3e5), fontproperties=font2, arrowprops=dict(facecolor=palette[2], shrink=0.05))
#annie.set_color(palette[2])
#plt.locator_params(axis='y',tight=True, numticks=4)
#plt.locator_params(axis='x',tight=True, nbins=4)
#plt.xticks(fontproperties=font2)
#plt.yticks(fontproperties=font2)
#main_ax.axes.xaxis.set_ticklabels([])
#main_ax.axes.yaxis.set_ticklabels([])
#main_ax.set_xlabel("Urine production", fontproperties=font2, labelpad=2)
#main_ax.set_ylabel("Bacteria in bladder", fontproperties=font2, labelpad=2)
#plt.xticks(fontproperties=font2)
#plt.yticks(fontproperties=font2)
#plt.tight_layout(rect=[-0.043,-0.044,1.033,1.036])
#plt.savefig('urine_rate_3mt_var.svg')

plt.show()