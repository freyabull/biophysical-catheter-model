""" Module to make graphs from multiple catheter simulations. """

import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import seaborn as sns
import os
from matplotlib.font_manager import FontProperties

# Graphics settings
plt.rcParams["figure.figsize"] = (3.2,3.2)
plt.rcParams['xtick.major.pad']='4'
plt.rcParams['ytick.major.pad']='4'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = '8'
font = FontProperties()
font.set_name('serif')
font.set_size(8)
font2 = FontProperties()
font2.set_name('serif')
font2.set_size(7)
colours = np.array([(68,120,33), (141,211,95), (255,179,128), (255,221,85), (179,179,179)])/255
palette = sns.color_palette(colours)
sns.set_palette(palette)
sns.set_style("ticks")
plt.rcParams['font.family'] = 'serif'

# Location of results files
ldur_folder = '../ExploringCatheter/results/urine_rates'
ldur_files = os.listdir(ldur_folder)

hdur_folder = '../ExploringCatheter/results/hd_urine_rates'
hdur_files = os.listdir(hdur_folder)

ldul_folder = '../ExploringCatheter/results/catheter_lengths'
ldul_files = os.listdir(ldul_folder)

hdul_folder = '../ExploringCatheter/results/hd_catheter_lengths'
hdul_files = os.listdir(hdul_folder)

#  Threshold value for blockage
blockage_thresh = 1e7

# Process data

# Low diffusion, urine rate

# Arrays for results
N = len(ldur_files)
ld_urine_rates = np.zeros(N)
ldur_time_to_block = np.zeros(N)
# Open files and read in data
for i in range(N):
    # Read file
    filename = os.path.join(ldur_folder, ldur_files[i])
    print(filename)
    df = pd.read_csv(filename, index_col=0, header=None, skiprows=3)
    info = pd.read_csv(filename,  nrows=1)
    if df.empty: break
    # Extract data
    ld_urine_rates[i] = info['urine rate']
    t_len = int(info['simulation length']/info['print interval']) # number of time steps
    dt = float(info['print interval']/3600) # print interval (hrs)
    maxy = np.zeros(t_len) # max value for inside density at each timestep
    for j in range(0, t_len):
        inside = np.array(df.iloc[4*j+2])
        maxy[j] = max(inside.astype(float))
    # Find times at which inside blocks
    index = np.searchsorted(maxy,blockage_thresh, side='right')
    ldur_time_to_block[i] = dt*index if index!=t_len else np.nan

# High diffusion, urine rate

# Arrays for results
N = len(hdur_files)
hd_urine_rates = np.zeros(N)
hdur_time_to_block = np.zeros(N)
# Open files and read in data
for i in range(N):
    # Read file
    filename = os.path.join(hdur_folder, hdur_files[i])
    print(filename)
    df = pd.read_csv(filename, index_col=0, header=None, skiprows=3)
    info = pd.read_csv(filename,  nrows=1)
    if df.empty: break
    # Extract data
    hd_urine_rates[i] = info['urine rate']
    t_len = int(info['simulation length']/info['print interval']) # number of time steps
    dt = float(info['print interval']/3600) # print interval (hrs)
    maxy = np.zeros(t_len) # max value for inside density at each timestep
    for j in range(0, t_len):
        inside = np.array(df.iloc[4*j+2])
        maxy[j] = max(inside.astype(float))
    # Find times at which inside blocks
    index = np.searchsorted(maxy,blockage_thresh, side='right')
    hdur_time_to_block[i] = dt*index if index!=t_len else np.nan

# Low diffusion, urethral length

# Arrays for results
N = len(ldul_files)
ld_lengths = np.zeros(N)
ldul_time_to_block = np.zeros(N)
# Open files and read in data
for i in range(N):
    # Read file
    filename = os.path.join(ldul_folder, ldul_files[i])
    print(filename)
    df = pd.read_csv(filename, index_col=0, header=None, skiprows=3)
    info = pd.read_csv(filename,  nrows=1)
    if df.empty: break
    # Extract data
    ld_lengths[i] = info['catheter length']
    t_len = int(info['simulation length']/info['print interval']) # number of time steps
    dt = float(info['print interval']/3600) # print interval (hrs)
    maxy = np.zeros(t_len) # max value for inside density at each timestep
    for j in range(0, t_len):
        inside = np.array(df.iloc[4*j+2])
        maxy[j] = max(inside.astype(float))
    # Find times at which inside blocks
    index = np.searchsorted(maxy,blockage_thresh, side='right')
    ldul_time_to_block[i] = dt*index if index!=t_len else np.nan

# High diffusion, urethral length

# Arrays for results
N = len(ldul_files)
hd_lengths = np.zeros(N)
hdul_time_to_block = np.zeros(N)
# Open files and read in data
for i in range(N):
    # Read file
    filename = os.path.join(hdul_folder, hdul_files[i])
    print(filename)
    df = pd.read_csv(filename, index_col=0, header=None, skiprows=3)
    info = pd.read_csv(filename,  nrows=1)
    if df.empty: break
    # Extract data
    hd_lengths[i] = info['catheter length']
    t_len = int(info['simulation length']/info['print interval']) # number of time steps
    dt = float(info['print interval']/3600) # print interval (hrs)
    maxy = np.zeros(t_len) # max value for inside density at each timestep
    for j in range(0, t_len):
        inside = np.array(df.iloc[4*j+2])
        maxy[j] = max(inside.astype(float))
    # Find times at which inside blocks
    index = np.searchsorted(maxy,blockage_thresh, side='right')
    hdul_time_to_block[i] = dt*index if index!=t_len else np.nan

# Set up figure/axes
gs_kw = dict(width_ratios=[1,1], height_ratios=[1,0.1,1])
big_fig, ax = plt.subplots(3,2, sharex='col', sharey=False, gridspec_kw=gs_kw)

ax[0,0].scatter(ld_urine_rates, ldur_time_to_block/24,s=4)
ax[0,0].set_yticks([168,169])
ax[2,0].scatter(hd_urine_rates, hdur_time_to_block/24,s=4)
ax[2,0].set_yticks([2.5,3])
ax[0,1].scatter(ld_lengths, ldul_time_to_block/24,s=4)
ax[0,1].set_yticks([400,800])
ax[2,1].scatter(hd_lengths, hdul_time_to_block/24,s=4)
ax[2,1].set_yticks([4,8])

ax[2,0].set_xlabel('Urine rate (mm$^3$s$^{-1}$)', labelpad=1.5)
ax[2,1].set_xlabel('Urethral length (mm)', labelpad=1.5)
sns.despine()
#ax[1,0].set_title('(a) Varying urine production \nrate for $D_S=10^{-8}$', fontproperties=font2, va='top')
ax[1,0].axis('off')
#ax[1,1].set_title('(b) Varying urethral length \nfor $D_S=10^{-8}$', fontproperties=font2, va='top')
ax[1,1].axis('off')

big_fig.add_subplot(111,frameon=False)
plt.tick_params(labelcolor='none',which='both',top=False,bottom=False,left=False,right=False)
plt.ylabel('Time to thick biofilm (days)', labelpad=3.5, fontproperties=font)
plt.text(0.2,0.525, '(a) Varying urine production \nrate for $D_S=10^{-8}$', fontproperties=font2, ha='center', va='top')
plt.text(0.775,0.525, '(b) Varying urethral length \nfor $D_S=10^{-8}$', fontproperties=font2, ha='center', va='top')
plt.text(0.2,-0.17, '(c) Varying urine production \nrate for $D_S=10^{-4}$', fontproperties=font2, ha='center', va='top')
plt.text(0.775,-0.17, '(d) Varying urethral length \nfor $D_S=10^{-4}$', fontproperties=font2, ha='center', va='top')

plt.tight_layout(rect=[-0.13,-0.165,1.075,1.04])
plt.savefig("diffusion_comparison.pdf")
#plt.show()
