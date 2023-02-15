""" Module to make graphs from multiple catheter simulations. """

import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import seaborn as sns
import os
from matplotlib.font_manager import FontProperties
from brokenaxes import brokenaxes
import matplotlib.ticker as ticker

# Graphics settings
#plt.rcParams["figure.figsize"] = (3.2,3.2)
plt.rcParams["figure.figsize"] = (1.55,1.55)
plt.rcParams['xtick.major.pad']='4'
plt.rcParams['ytick.major.pad']='4'
plt.rcParams['text.usetex'] = True
font = FontProperties()
font.set_name('serif')
font.set_size(8)
colours = np.array([(68,120,33), (141,211,95), (255,179,128), (255,221,85), (179,179,179), (199,143,104), (144,106,81)])/255
palette = sns.color_palette(colours)
palette3 = sns.dark_palette(colours[2], n_colors=5,reverse=True)
sns.set_palette(palette)
sns.set_style("ticks")
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = '8'
print(palette[5])

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
fig_names = ['ic_bag_big.pdf' ,'ic_bag.pdf', 'ic_bag_3.pdf', 'ic_bladder_big.pdf', 'ic_bladder.pdf', 'ic_bladder_3.pdf', 'ic_skin_big.pdf', 'ic_skin.pdf', 'ic_skin_3.pdf', 'ic_uniform_big.pdf', 'ic_uniform.pdf', 'ic_uniform_3.pdf']
fig_svgs = ['ic_bag.svg', 'ic_bladder.svg', 'ic_skin.svg', 'ic_uniform.svg']
figs = []
for file in files:
    df = pd.read_csv(file, index_col=0, header=None, skiprows=3)
    info = pd.read_csv(file,  nrows=1)
    t_len = int(info['simulation length']/info['print interval']) # number of time steps
    dt = float(info['print interval']/3600) # print interval (hrs)
    time = dt*np.array(range(0, t_len)) # time series
    catheter_length = float(info['catheter length'])
    x_len = int(info['print num steps']) # number of x steps
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

    plt.rcParams["figure.figsize"] = (3.2,3.2)
    fig, ax = plt.subplots()
    plt.tight_layout(rect=[-0.003,0.005,1.055,1.02])
    plt.xticks(fontproperties=font)
    ax.spines['left'].set_position('zero')
    maxx = max(max(outside[-1]),max(inside[-1]))
    bax = brokenaxes(xlims=((-1.03*maxx,0),(0,1.03*maxx)))
    ax.axvspan(0,1,color=palette[4])
    inside_50 = inside[index24*50][::-1]
    outside_50 = outside[index24*50]
    inside_60 = inside[index24*60][::-1]
    outside_60 = outside[index24*60]
    bax.plot(inside_50,x,color=palette[1],linestyle='dashed',linewidth=1)
    bax.plot(-1*outside_50,x,color=palette[1],linestyle='dashed',linewidth=1)
    bax.plot(inside_60,x,color=palette[0],linewidth=2)
    bax.plot(-1*outside_60,x,color=palette[0],linewidth=2)
    ax.axis('off')
    bax.first_col[0].yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
    bax.set_ylabel('Distance up urethra (mm)', labelpad=20)
    bax.set_xlabel('Bacterial density (mm$^{-2}$)', labelpad=15)
    #bax.last_row[0].set_xticklabels(['','1','','0'])
    #bax.last_row[1].set_xticklabels(['','0','','1',''])
    bax.last_row[0].set_xticklabels(['','$10^7$','','0'])
    bax.last_row[1].set_xticklabels(['','0','','$10^7$',''])
    ax.text(0.055,1.001,'Extraluminal')
    ax.text(0.665,1.001, 'Intraluminal')
    figs.append(fig)
    plt.savefig(fig_names.pop(0))
    #plt.savefig(fig_svgs.pop(0))



    plt.rcParams["figure.figsize"] = (1.55,1.55)
    fig, ax = plt.subplots()
    plt.tight_layout(rect=[0.005,0.01,1.12,1.035])
    plt.xticks(fontproperties=font)
    ax.spines['left'].set_position('zero')
    maxx = max(max(outside[-1]),max(inside[-1]))
    bax = brokenaxes(xlims=((-1.03*maxx,0),(0,1.03*maxx)))
    #bax.axvspan(0,maxx,color=palette[3], zorder=1, alpha=0.35)
    ax.axvspan(0,1,color=palette[4])
    #bax.fill_betweenx(x,inside[index24*180][::-1],color=palette[0], zorder=2, label='{} days'.format(8))
    #bax.fill_betweenx(x,-1*outside[index24*180],color=palette[0], zorder=2)
    #bax.fill_betweenx(x,inside[index24*60][::-1],color=palette[1], zorder=3, label='{} days'.format(3))
    #bax.fill_betweenx(x,-1*outside[index24*60],color=palette[1], zorder=3)
    #bax.fill_betweenx(x,inside[index24*15][::-1],color=palette[2], zorder=4, label='{} days'.format(1))
    #bax.fill_betweenx(x,-1*outside[index24*15],color=palette[2], zorder=4)
    inside_50 = inside[index24*50][::-1]
    #inside_50[ inside_50<1e-30 ] = np.nan
    outside_50 = outside[index24*50]
    inside_60 = inside[index24*60][::-1]
    outside_60 = outside[index24*60]
    #outside_60[ outside_60<1e-30 ] = np.nan
    bax.plot(inside_50,x,color=palette[1],linestyle='dashed',linewidth=1)
    bax.plot(-1*outside_50,x,color=palette[1],linestyle='dashed',linewidth=1)
    bax.plot(inside_60,x,color=palette[0],linewidth=2)
    bax.plot(-1*outside_60,x,color=palette[0],linewidth=2)
    ax.axis('off')
    #bax.first_col[0].yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
    #bax.set_ylabel('Distance up urethra (mm)', labelpad=20,position=(0,0.43))
    bax.set_ylabel('$\it{x}$ (mm)', labelpad=20)
    bax.set_xlabel('Bacterial density (mm$^{-2}$)', labelpad=15,position=(0.40,0))
    #bax.last_row[0].set_xticklabels(['','1','','0'])
    #bax.last_row[1].set_xticklabels(['','0','','1',''])
    bax.last_row[0].set_xticklabels(['','$10^7$',''])
    bax.last_row[1].set_xticklabels(['','0','$10^7$',''])
    #bax.text(0.19,0.9,r'$\times 10^7$')
    #bax.text(-1.2,45,'Outside')
    #bax.text(0.2,45,'Inside')
    ax.text(-0.012,1.01,'Outside')
    ax.text(0.635,1.01, 'Inside')
    figs.append(fig)
    plt.savefig(fig_names.pop(0))
    plt.savefig(fig_svgs.pop(0))

    plt.rcParams["figure.figsize"] = (1.55,1.55)
    fig, (ax1, ax2) = plt.subplots(2,1, sharex=True, sharey=True)
    ax2.plot(x,inside_50,color=palette[1],linestyle='dashed',linewidth=1)
    ax1.plot(x,outside_50,color=palette[1],linestyle='dashed',linewidth=1)
    ax2.plot(x,inside_60,color=palette[0],linewidth=2)
    ax1.plot(x,outside_60,color=palette[0],linewidth=2)
    sns.despine()
    ax1.set_title('Extraluminal', fontproperties=font, pad=3.0)
    ax2.set_title('Intraluminal', fontproperties=font, pad=3.0)
    ax1.set_ylim(0,1.035*maxx)
    ax2.set_ylim(0,1.035*maxx)
    ax1.set_yticklabels(['0', '$10^7$'])
    ax2.set_yticklabels(['0', '$10^7$'])
    fig.add_subplot(111,frameon=False)
    plt.tick_params(labelcolor='none',which='both',top=False,bottom=False,left=False,right=False)
    plt.ylabel('Bacterial density (mm$^{-2}$)', labelpad=1.5, position=(0,0.47))
    plt.xlabel('Distance up urethra (mm)', labelpad=1.5, position=(0.3,0))
    plt.tight_layout(rect=[-0.24,-0.195,1.13,1.12])
    figs.append(fig)
    plt.savefig(fig_names.pop(0))

#plt.show()

# One giant figure attempt
plt.rcParams["figure.figsize"] = (3.2,3.4)
plt.rcParams['xtick.major.pad']='4'
plt.rcParams['ytick.major.pad']='4'
plt.rcParams['font.size'] = 8.0
print(plt.rcParams['font.family'])
font = FontProperties()
font.set_name('serif')
font.set_size(8)
font2 = FontProperties()
font2.set_name('serif')
font2.set_size(7)
colours = np.array([(68,120,33), (141,211,95), (255,179,128), (255,221,85), (179,179,179)])/255
#palette = sns.color_palette(colours)
#sns.set_palette(palette)
sns.set_style("ticks")
plt.rcParams['font.family'] = 'serif'
gs_kw = dict(width_ratios=[1,0.05,1], height_ratios=[1,1,0.1,1,1])
big_fig, ax = plt.subplots(5,3, sharex=True, sharey=True, gridspec_kw=gs_kw)

bigfiles = ['../BasicCatheter/initial_conditions/skin.csv','../BasicCatheter/initial_conditions/bag.csv','../BasicCatheter/initial_conditions/uniform.csv','../BasicCatheter/initial_conditions/bladder.csv']

for f in range(4):
    file = bigfiles[f]
    df = pd.read_csv(file, index_col=0, header=None, skiprows=3)
    info = pd.read_csv(file,  nrows=1)
    t_len = int(info['simulation length']/info['print interval']) # number of time steps
    dt = float(info['print interval']/3600) # print interval (hrs)
    time = dt*np.array(range(0, t_len)) # time series
    catheter_length = float(info['catheter length'])
    x_len = int(info['print num steps']) # number of x steps
    dx = catheter_length/(x_len-1) # x step interval (mm)
    x = dx*np.array(range(0,x_len)) # x series
    outside = np.zeros((t_len,df.shape[1]))
    inside = np.zeros((t_len,df.shape[1]))
    for i in range(0, t_len):
        outside[i] = df.iloc[4*i]
        inside[i] = df.iloc[4*i+2]
    outside = outside*1e-7
    inside = inside*1e-7
    maxx = max(max(outside[-1]),max(inside[-1]))
    index24 = int(24/dt) #find index for 24 hrs
    inside_50 = inside[index24*50][::-1]
    outside_50 = outside[index24*50]
    inside_60 = inside[index24*60][::-1]
    outside_60 = outside[index24*60]
    #print(plt.rcParams['font.family'])

    i = (f//2)*3
    j = f%2*2
    print(i,j)
    ax[i+1,j].plot(x,inside_50,color=palette[2],linestyle='dashed',linewidth=1)
    ax[i,j].plot(x,outside_50,color=palette[1],linestyle='dashed',linewidth=1)
    ax[i+1,j].plot(x,inside_60,color=palette3[2],linewidth=2)
    ax[i,j].plot(x,outside_60,color=palette[0],linewidth=2)
    sns.despine()
    #ax[i,j].set_title('Extraluminal', fontproperties=font, pad=3.0)
    #ax[i+1,j].set_title('Intraluminal', fontproperties=font, pad=3.0)
    ax[i,j].set_ylim(0,1.035*maxx)
    ax[i+1,j].set_ylim(0,1.035*maxx)
    ax[i,j].set_yticks([0,1])
    ax[i+1,j].set_yticks([0,1])
    ax[i,j].set_yticklabels(['0', '$10^7$'])
    ax[i+1,j].set_yticklabels(['0', '$10^7$'])
    #print(plt.rcParams['font.family'])


#ax[2,0].set_title('(a) Infection originates \nfrom skin', fontproperties=font2, va='top')
#ax[2,0].set_title('(a) Skin', fontproperties=font)
ax[2,0].axis('off')
ax[2,1].axis('off')
ax[2,2].axis('off')
ax[0,1].axis('off')
ax[1,1].axis('off')
ax[3,1].axis('off')
ax[4,1].axis('off')
#ax[2,1].set_title('(b) Infection originates \nfrom drainage bag', fontproperties=font2, va='top')
#ax[2,1].set_title('(b) Drainage bag', fontproperties=font)

big_fig.add_subplot(111,frameon=False)
plt.tick_params(labelcolor='none',which='both',top=False,bottom=False,left=False,right=False)
plt.ylabel('Bacterial density (mm$^{-2}$)', labelpad=3.0, fontproperties=font)
plt.xlabel('Distance up urethra (mm)', labelpad=1.5, fontproperties=font)
#plt.text(0.225,-0.13, '(c) Initial contamination \nis uniform over \nextraluminal surface', fontproperties=font2, ha='center', va='top')
#plt.text(0.775,-0.13, '(d) Pre-existing bladder \ncontamination before \ncatheter insertion', fontproperties=font2, ha='center', va='top')
plt.text(0.0,1.02, '(a) Skin', fontproperties=font)
plt.text(0.575,1.02, '(b) Drainage bag', fontproperties=font)
plt.text(0.0,0.45, '(c) Insertion', fontproperties=font)
plt.text(0.575,0.45, '(d) Bladder', fontproperties=font)
#plt.text(1.04,0.9, 'Ext.', fontproperties=font, rotation='vertical', va='center')
#plt.text(1.04,0.65, 'Int.', fontproperties=font, rotation='vertical', va='center')
#plt.text(1.04,0.33, 'Ext.', fontproperties=font, rotation='vertical', va='center')
#plt.text(1.04,0.08, 'Int.', fontproperties=font, rotation='vertical', va='center')
plt.tight_layout(rect=[-0.12,-0.093,1.05,1.03])
plt.savefig("ic_bigfig.pdf")



