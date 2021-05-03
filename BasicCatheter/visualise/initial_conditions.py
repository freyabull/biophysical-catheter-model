""" Module to make graphs from multiple catheter simulations. """

import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import seaborn as sns
import os
from matplotlib.font_manager import FontProperties
from brokenaxes import brokenaxes

# Graphics settings
#plt.rcParams["figure.figsize"] = (11/2.54,8/2.54)
plt.rcParams["figure.figsize"] = (3.2,3.2)
plt.rcParams['xtick.major.pad']='4'
plt.rcParams['ytick.major.pad']='4'
font = FontProperties()
font.set_name('serif')
font.set_size(8)
print(plt.rcParams.keys())
plt.rcParams['font.size'] = 8.0
colours = np.array([(68,120,33), (141,211,95), (255,179,128), (255,221,85), (179,179,179)])/255
palette = sns.color_palette(colours)
sns.set_palette(palette)
sns.set_style("ticks")
plt.rcParams['font.family'] = 'serif'

## Location of results files
#folder = '../BasicCatheter/initial_conditions'

## Actually lets open files one at a time
## Open one file and manipulate it as a demo then repeat for other cases. Create 4 medium figures
## Open file
#file = '../BasicCatheter/initial_conditions/skin.csv'
#df = pd.read_csv(file, index_col=0, header=None, skiprows=3)
#info = pd.read_csv(file,  nrows=1)

## Read outside and inside
##Extract data series from dataframe
#t_len = int(info['simulation length']/info['print interval']) # number of time steps
#dt = float(info['print interval']/3600) # print interval (hrs)
#time = dt*np.array(range(0, t_len)) # time series
#catheter_length = float(info['catheter length'])
#x_len = int(info['num of x steps']) # number of x steps
#dx = catheter_length/(x_len-1) # x step interval (mm)
#x = dx*np.array(range(0,x_len)) # x series
#outside = np.zeros((t_len,df.shape[1]))
#inside = np.zeros((t_len,df.shape[1]))
#for i in range(0, t_len):
#    outside[i] = df.iloc[4*i]
#    inside[i] = df.iloc[4*i+2]

##find index for 24 hrs
#index24 = int(24/dt)
## Plot distributions at 3 time points

#fig3, axy = plt.subplots()
#axy.spines['left'].set_position('zero')
#bax = brokenaxes(xlims=((-1*max(outside[-1]),0),(0,max(inside[-1]))))
##bax.plot(inside[-1][::-1],x,color=palette[0])
##bax.plot(-1*outside[-1],x,color=palette[0])
##bax.plot(inside[index24*5][::-1],x,color=palette[1])
##bax.plot(-1*outside[index24*5],x,color=palette[1])
##bax.plot(inside[index24*2][::-1],x,color=palette[2])
##bax.plot(-1*outside[index24*2],x,color=palette[2])
#bax.axvspan(0,max(inside[-1]),color=palette[3], zorder=1)
#axy.axvspan(0,1,color=palette[4])
#bax.fill_betweenx(x,inside[-1][::-1],color=palette[0], zorder=2, label='{} days'.format(int(info['simulation length']/(3600*24))))
#bax.fill_betweenx(x,-1*outside[-1],color=palette[0], zorder=2)
#bax.fill_betweenx(x,inside[index24*4][::-1],color=palette[1], zorder=3, label='{} days'.format(4))
#bax.fill_betweenx(x,-1*outside[index24*4],color=palette[1], zorder=3)
#bax.fill_betweenx(x,inside[index24*2][::-1],color=palette[2], zorder=4, label='{} days'.format(2))
#bax.fill_betweenx(x,-1*outside[index24*2],color=palette[2], zorder=4)
#axy.axis('off')
#bax.set_ylabel('Distance up urethra (mm)', labelpad=22)
#bax.set_xlabel('Bacterial density (mm$^{-3}$)')
#bax.legend()


#plt.show()


# Full attempt
files = ['../BasicCatheter/initial_conditions/bag.csv','../BasicCatheter/initial_conditions/bladder.csv','../BasicCatheter/initial_conditions/skin.csv','../BasicCatheter/initial_conditions/uniform.csv']
fig_names = ['ic_bag.pdf', 'ic_bladder.pdf', 'ic_skin.pdf', 'ic_uniform.pdf']
fig_svgs = ['ic_bag.svg', 'ic_bladder.svg', 'ic_skin.svg', 'ic_uniform.svg']
figs = []
for file in files:
    df = pd.read_csv(file, index_col=0, header=None, skiprows=3)
    info = pd.read_csv(file,  nrows=1)
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
    outside = outside*1e-7
    inside = inside*1e-7
    index24 = int(24/dt) #find index for 24 hrs
    fig, ax = plt.subplots()
    plt.tight_layout(rect=[0.02,0.005,1.06,1.045])
    plt.xticks(fontproperties=font)
    ax.spines['left'].set_position('zero')
    bax = brokenaxes(xlims=((-1*max(outside[-1]),0),(0,max(inside[-1]))))
    bax.axvspan(0,max(inside[-1]),color=palette[3], zorder=1, alpha=0.35)
    ax.axvspan(0,1,color=palette[4])
    bax.fill_betweenx(x,inside[index24*8][::-1],color=palette[0], zorder=2, label='{} days'.format(8))
    bax.fill_betweenx(x,-1*outside[-1],color=palette[0], zorder=2)
    bax.fill_betweenx(x,inside[index24*4][::-1],color=palette[1], zorder=3, label='{} days'.format(4))
    bax.fill_betweenx(x,-1*outside[index24*4],color=palette[1], zorder=3)
    bax.fill_betweenx(x,inside[index24*2][::-1],color=palette[2], zorder=4, label='{} days'.format(2))
    bax.fill_betweenx(x,-1*outside[index24*2],color=palette[2], zorder=4)
    ax.axis('off')
    bax.set_ylabel('Distance up urethra (mm)', labelpad=25)
    bax.set_xlabel('Bacterial density (mm$^{-3}$)', labelpad=15)
    #bax.legend()
    bax.text(5.8,-5.0,r'$\times 10^7$')
    figs.append(fig)
    #plt.savefig(fig_names.pop(0))
    plt.savefig(fig_svgs.pop(0))

#plt.show()
