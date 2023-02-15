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
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = '8'
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
folder = '../ExploringCatheter/results/sump_volumes'
files = os.listdir(folder)

# Arrays for results
N = len(files)
sump_volumes = np.zeros(N)
times = np.zeros(N)
max_infect = np.zeros(N)
time_to_max = np.zeros(N)

#  Threshold value for blockage
thresh = 1e7


# Open files and read in data
for i in range(N):
    filename = os.path.join(folder, files[i])
    print(filename)
    # Read data
    df = pd.read_csv(filename, index_col=0, header=None, skiprows=3)
    info = pd.read_csv(filename,  nrows=1)
    if df.empty: break

    growth_rate = float(info['bladder growth rate'])
    urine_rate = float(info['urine rate'])
    sump_volumes[i] = info['sump volume']
    t_len = int(info['simulation length']/info['print interval']) # number of time steps
    dt = float(info['print interval']/3600) # print interval (hrs)

    bladder = np.zeros(t_len)
    inside = np.zeros((t_len,df.shape[1]))
    maxy = np.zeros(t_len)
    for j in range(0, t_len):
        bladder[j] = df.iloc[4*j+1][1]
        inside[j] = df.iloc[4*j+2]
        maxy[j] = max(inside[j])

    index = np.searchsorted(maxy,thresh, side='right')
    times[i] = dt*index if index!=t_len else np.nan

    # Find maximum bladder infection & time to infection
    max_infect[i] = np.max(bladder)
    time_to_max[i] = dt*np.searchsorted(bladder,0.99*max_infect[i], side='right')

plt.scatter(sump_volumes[sump_volumes<=1e6],times[sump_volumes<=1e6],s=4, color=palette[0])
sns.despine()
plt.xscale('log')
plt.ylim(4033,4073)
plt.locator_params(axis='x',tight=True, numticks=4)
plt.locator_params(axis='y',tight=True, nbins=3)
plt.xticks(fontproperties=font)
plt.yticks(fontproperties=font)
plt.xlabel('Urine volume (mm$^3$)', fontproperties=font, labelpad=1.5, position=(0.46,0))
plt.ylabel('Time to thick biofilm (hr)', fontproperties=font, labelpad=2,position=(0,0.41))
plt.tight_layout(rect=[-0.072,-0.112,1.08,1.07])
plt.savefig('sump_volumes_blockages.pdf')
#plt.show()

plt.rcParams["figure.figsize"] = (1.8,2.5)
plt.rcParams['text.usetex'] = True
plt.rcParams['xtick.major.pad']='4'
plt.rcParams['ytick.major.pad']='4'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = '12'
font = FontProperties()
font.set_name('serif')
font.set_size(12)


#plt.figure()
#plt.scatter(sump_volumes[sump_volumes<=1e6],times[sump_volumes<=1e6],s=4, color=palette[0])
#sns.despine()
#plt.xscale('log')
#plt.ylim(4033,4073)
#plt.locator_params(axis='x',tight=True, numticks=4)
#plt.locator_params(axis='y',tight=True, nbins=3)
#plt.xticks(fontproperties=font)
#plt.yticks(fontproperties=font)
#plt.xlabel('Urine volume (mm$^3$)', fontproperties=font, labelpad=1.5)#, position=(0.46,0))
#plt.ylabel('Time to thick biofilm (hr)', fontproperties=font, labelpad=2)#,position=(0,0.41))
#plt.tight_layout(rect=[-0.05,-0.05,1.08,1.07])
#plt.savefig('thesis_sump_volumes_d-8.pdf')
#plt.show()

#plt.rcParams["figure.figsize"] = (3.6,2.5)
plt.figure()
#fig, (ax1, ax2) = plt.subplots(1,2, sharey=True)
ax1 = plt.subplot()
print('sump volumes ', sump_volumes) # how big is 1L? this should be absolute max sump volume, maybe even 500 mL
print('max infect ', max_infect)
ax1.scatter(sump_volumes[sump_volumes>=1e4]*1e-5,max_infect[sump_volumes>=1e4]*1e-5,zorder=2,color=palette[0])
ax1.axvline(urine_rate*1e-5/growth_rate, color=palette[1], linestyle='--')
sns.despine()
plt.locator_params(tight=True, nbins=3)
ax1.set_xlabel('$V$ ($10^5$ mm$^3$)', labelpad=1.5)
ax1.set_ylabel('Density ($10^5$ mm$^{-3}$)', labelpad=2)
plt.tight_layout(rect=[-0.09,-0.07,1.09,1.07])
plt.savefig('thesis_sump_volumes.pdf')

plt.figure()
ax2 = plt.subplot()
ax2.scatter(urine_rate*1e4/sump_volumes[sump_volumes>=1e4],max_infect[sump_volumes>=1e4]*1e-5,zorder=2,color=palette[0])
ax2.axvline(growth_rate*1e4, color=palette[1], linestyle='--')
sns.despine()
plt.xticks(fontproperties=font)
ax2.set_yticklabels([])
plt.locator_params(tight=True, nbins=3)
ax2.set_xlabel('$\lambda/V$ ($10^{-4}$ s$^{-1}$)', labelpad=1.5)
#plt.ylabel('Density ($10^5$ mm$^{-3}$)', labelpad=2)
#plt.tight_layout(rect=[-0.04,-0.06,1.06,1.065], pad=1.0)
plt.tight_layout(rect=[-0.09,-0.07,1.09,1.05])
plt.savefig('thesis_sump_volumes_i.pdf')
#plt.savefig('thesis_sump_volumes_2.pdf')

#plt.rcParams["figure.figsize"] = (1.8,2.5)
plt.figure()
ax2 = plt.subplot()
ax2.scatter(sump_volumes[sump_volumes>=1e4]*1e-5,time_to_max[sump_volumes>=1e4]/24,zorder=2,color=palette[0])
plt.axvline(urine_rate*1e-5/growth_rate, color=palette[1], linestyle='--')
sns.despine()
plt.locator_params(tight=True, nbins=3)
plt.xlabel('$V$ ($10^5$ mm$^3$)', labelpad=1.5)
plt.ylabel('Time (days)', labelpad=2)
plt.tight_layout(rect=[-0.09,-0.07,1.09,1.06])
plt.savefig('thesis_sump_volumes_t.pdf')
plt.show()


