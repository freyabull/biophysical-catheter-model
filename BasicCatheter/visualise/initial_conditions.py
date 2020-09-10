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
folder = '../BasicCatheter/initial_conditions'

# Actually lets open files one at a time
# Open one file and manipulate it as a demo then repeat for other cases. Create 4 medium figures
# Open file
file = '../BasicCatheter/initial_conditions/skin.csv'
df = pd.read_csv(file, index_col=0, header=None, skiprows=3)
info = pd.read_csv(file,  nrows=1)

# Read outside and inside
#Extract data series from dataframe
t_len = int(info['simulation length']/info['print interval']) # number of time steps
dt = float(info['print interval']/3600) # print interval (hrs)
time = dt*np.array(range(0, t_len)) # time series
catheter_length = float(info['catheter length'])
x_len = int(info['num of x steps']) # number of x steps
dx = catheter_length/(x_len-1) # x step interval (mm)
x = dx*np.array(range(0,x_len)) # x series
outside = np.zeros((t_len,df.shape[1]))
inside = np.zeros((t_len,df.shape[1]))
for i in range(0, t_len):
    outside[i] = df.iloc[4*i]
    inside[i] = df.iloc[4*i+2]

# Plot distributions at 3 time points
fig = plt.figure(constrained_layout=True)
gs = fig.add_gridspec(1, 3)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1:])
ax1.plot(outside[-1],x)
ax1.set_xlim(1e7,0)
ax2.plot(inside[-1],x)
ax1.invert_yaxis() 
ax2.set_xlim(0,2e7)
sns.despine(ax=ax2)
sns.despine(ax=ax1, left=True, right = False)
ax2.set_ylabel('Distance from catheter tip', labelpad=2)
ax1.yaxis.tick_right()
print(ax1.get_ylim())
print(ax1.get_yticks())
print(ax2.get_ylim())
print(ax2.get_yticks())
#ax2.set_ylim(ax1.get_ylim()[1],ax1.get_ylim()[0])
#ax2.set_yticks(catheter_length-ax1.get_yticks())
ax2.axes.yaxis.set_ticklabels([])
print(ax1.get_ylim())
print(ax1.get_yticks())
print(ax2.get_ylim())
print(ax2.get_yticks())

# Lets start again. CAn I plot both on the same axis (one +ve, one -ve)?


plt.show()