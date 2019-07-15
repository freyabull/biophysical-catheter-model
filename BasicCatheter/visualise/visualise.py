""" Module to make graphs to visualise the catheter simulation. """

import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import seaborn as sns
from textwrap import wrap


# Read data into pandas dataframe
df = pd.read_csv('../BasicCatheter/results.csv', index_col=0, header=None, skiprows=3)
info = pd.read_csv('../BasicCatheter/results.csv',  nrows=1)


#Extract data series from dataframe
t_len = int(info['simulation length']/info['print step']) # number of time steps
dt = float(info['print step']/3600) # print interval
time = dt*np.array(range(0, t_len)) # time series

x_len = int(info['num of x steps']) # number of x steps
dx = float(info['catheter length']/(x_len-1)) # x step interval
x = dx*np.array(range(0,x_len)) # x series

outside = np.zeros((t_len,df.shape[1]))
bladder = np.zeros(t_len)
inside = np.zeros((t_len,df.shape[1]))
for i in range(0, t_len):
    outside[i] = df.iloc[3*i]
    bladder[i] = df.iloc[3*i+1][1]
    inside[i] = df.iloc[3*i+2]


# Data processing required for multiple graphs
# Find locations of peak concentration
peak_outside = np.argmax(outside, axis=1)
peak_inside = np.argmax(inside, axis=1)

# Find concentration values at peaks
opeaks = np.array([outside[i][peak_outside[i]] for i in range(1,t_len)])
ipeaks = np.array([inside[i][peak_inside[i]] for i in range(1,t_len)])


# Set up figure/environment
plt.rc('text', usetex=True) # Enable latex
sns.set_palette(sns.dark_palette("sky blue", n_colors=8, input="xkcd", reverse=True)) # Set colour palette
fig, ((ax_bt, ax_op, ax_ipe, ax_ipl), (ax_ow, ax_iw, ax_owf, ax_iwf)) = plt.subplots(2,4) # Create figure with subplots
fig.suptitle('\n'.join(wrap('Model parameters: {}'.format([(i,info[i][0]) for i in info]),300)), fontsize=10) # Add model parameters to figure


# Figure 1
# Bladder concentration against time
min_bladder = float(1.0/(1000.0*info['sump volume'])) # Find reasonable minimum bladder concentration
bladder = bladder*(bladder>min_bladder) # Remove unphysical small values
ax_bt.plot(time, bladder)
ax_bt.xaxis.set_major_locator(plt.MultipleLocator(12)) # Format time axis ticks
ax_bt.set_yscale('log')
ax_bt.set_xlabel('Time (hrs)')
ax_bt.set_ylabel('Bacterial concentration (mm$^{-3}$)')
ax_bt.set_title('Bladder')


# Figure 2
# Superposition of extraluminal wave profiles
otimes = np.unique([np.argmax(opeaks>=opeaks[-1]*0.2*i) for i in range(1,5)]) # Find interesting times to plot: peak concentration = (0.2,0.4,0.6,0.8)*max
for i in otimes: ax_op.plot(x,outside[i], label='$t = {}$ hrs'.format(time[i])) # Plot wave profile for times found above
ax_op.plot(x,outside[-1], label='$t = {}$ hrs'.format(time[-1])) # Plot wave profile for final time (ie peak concentration = max)
ax_op.ticklabel_format(axis='y', style='sci', scilimits=(0,0)) # Axis labels for concentration in standard form
ax_op.set_xlabel('Distance from catheter base (mm)')
ax_op.set_ylabel('Bacterial concentration (mm$^{-2}$)')
ax_op.set_title('Extraluminal wave profiles')
ax_op.legend(loc='upper right')


# Figure 3
# Superposition of intraluminal wave profiles at early times
istimes = np.unique([np.argmax(ipeaks>=(20*i)) for i in range(1,6)]) # Find interesting times to plot: concentration = (20,40,60,80,100)
for i in istimes: ax_ipe.plot(x,inside[i], label='$t = {}$ hrs'.format(time[i])) # Plot wave profile for times found above
ax_ipe.set_xlabel('Distance from catheter tip (mm)')
ax_ipe.set_ylabel('Bacterial concentration (mm$^{-2}$)')
ax_ipe.set_title('Intraluminal wave profiles, early time')
ax_ipe.legend(loc='upper right')


# Figure 4
# Superposition of intraluminal wave profiles at late times
itimes = np.unique([np.argmax(ipeaks>=ipeaks[-1]*0.2*i) for i in range(1,5)]) # Find interesting times to plot: peak concentration = (0.2,0.4,0.6,0.8)*max
for i in itimes: ax_ipl.plot(x,inside[i], label='$t = {}$ hrs'.format(time[i])) # Plot wave profile for times found above
ax_ipl.plot(x,inside[-1],label='$t = {}$ hrs'.format(time[-1])) # Plot wave profile for final time (ie peak concentration = max)
ax_ipl.ticklabel_format(axis='y', style='sci', scilimits=(0,0)) # Axis labels for concentration in standard form
ax_ipl.set_xlabel('Distance from catheter tip (mm)')
ax_ipl.set_ylabel('Bacterial concentration (mm$^{-2}$)')
ax_ipl.set_title('Intraluminal wave profiles, late time')
ax_ipl.legend(loc='upper right')


# Figure 5
# Location of extraluminal peak concentration with measure of peak width
# Find width at 0.99*peak
po_min = np.empty(t_len)
po_max = np.empty(t_len)
po_min[:] = np.nan
po_max[:] = np.nan
for i in range(1,t_len):
    po_range = np.nonzero(outside[i]>=(outside[i][int(peak_outside[i])]*0.99)) # Find points within peak
    if len(po_range[0]): po_min[i] = po_range[0][0]; po_max[i] = po_range[0][-1] # Find edges of peak
    pass
ax_ow.plot(time[1:], peak_outside[1:]*dx) # Plot location of peak concentration
ax_ow.fill_between(time, po_min*dx, po_max*dx, alpha=0.2) # Add shading indicating peak width
ax_ow.xaxis.set_major_locator(plt.MultipleLocator(12)) # Format time axis ticks
ax_ow.set_xlabel('Time (hrs)')
ax_ow.set_ylabel('Distance from catheter base (mm)')
ax_ow.set_title('Position of extraluminal peak concentration')


# Figure 6
# Location of intraluminal peak concentration with measure of peak width
# Find width at 0.99*peak
pi_min = np.empty(t_len)
pi_max = np.empty(t_len)
pi_min[:] = np.nan
pi_max[:] = np.nan
for i in range(1,t_len):
    pi_range = np.nonzero(inside[i]>=(inside[i][int(peak_inside[i])]*0.99)) # Find points within peak
    if len(pi_range[0]): pi_min[i] = pi_range[0][0]; pi_max[i] = pi_range[0][-1] # Find edges of peak
    pass
ax_iw.plot(time[1:], peak_inside[1:]*dx) # Plot location of peak concentration
ax_iw.fill_between(time, pi_min*dx, pi_max*dx, alpha=0.2) # Add shading indicating peak width
ax_iw.xaxis.set_major_locator(plt.MultipleLocator(12)) # Format time axis ticks
ax_iw.set_xlabel('Time (hrs)')
ax_iw.set_ylabel('Distance from catheter tip (mm)')
ax_iw.set_title('Position of intraluminal peak concentration')


# Figure 7
# Extraluminal wavefronts at fixed concentration
# Find locations of wavefronts
ofrontl = np.empty([6,t_len])
ofrontu = np.empty([6,t_len])
ofrontl[:][:] = np.nan
ofrontu[:][:] = np.nan
for i in range(1,t_len):
    ofront0 = np.nonzero(outside[i]>=1e0)
    if len(ofront0[0]): ofrontl[0][i] = ofront0[0][0]; ofrontu[0][i] = ofront0[0][-1]
    ofront1 = np.nonzero(outside[i]>=1e1)
    if len(ofront1[0]): ofrontl[1][i] = ofront1[0][0]; ofrontu[1][i] = ofront1[0][-1]
    ofront2 = np.nonzero(outside[i]>=1e2)
    if len(ofront2[0]): ofrontl[2][i] = ofront2[0][0]; ofrontu[2][i] = ofront2[0][-1]
    ofront3 = np.nonzero(outside[i]>=1e3)
    if len(ofront3[0]): ofrontl[3][i] = ofront3[0][0]; ofrontu[3][i] = ofront3[0][-1]
    ofront4 = np.nonzero(outside[i]>=1e4)
    if len(ofront4[0]): ofrontl[4][i] = ofront4[0][0]; ofrontu[4][i] = ofront4[0][-1]
    ofront5 = np.nonzero(outside[i]>=1e5)
    if len(ofront5[0]): ofrontl[5][i] = ofront5[0][0]; ofrontu[5][i] = ofront5[0][-1]
    pass
for i in range(6): ax_owf.fill_between(time, ofrontl[i]*dx, ofrontu[i]*dx, label='$n = 10^{}$ mm$^{{-2}}$'.format(i))
ax_owf.xaxis.set_major_locator(plt.MultipleLocator(12)) # Format time axis ticks
ax_owf.set_xlabel('Time (hrs)')
ax_owf.set_ylabel('Distance from catheter base (mm)')
ax_owf.set_title('Wavefront outside')
ax_owf.legend(loc='upper right')


# Figure 8
# Intraluminal wavefronts at fixed concentration
# Find locations of wavefronts
ifrontl = np.empty([6,t_len])
ifrontu = np.empty([6,t_len])
ifrontl[:][:] = np.nan
ifrontu[:][:] = np.nan
for i in range(1,t_len):
    ifront0 = np.nonzero(inside[i]>=1e0)
    if len(ifront0[0]): ifrontl[0][i] = ifront0[0][0]; ifrontu[0][i] = ifront0[0][-1]
    ifront1 = np.nonzero(inside[i]>=1e1)
    if len(ifront1[0]): ifrontl[1][i] = ifront1[0][0]; ifrontu[1][i] = ifront1[0][-1]
    ifront2 = np.nonzero(inside[i]>=1e2)
    if len(ifront2[0]): ifrontl[2][i] = ifront2[0][0]; ifrontu[2][i] = ifront2[0][-1]
    ifront3 = np.nonzero(inside[i]>=1e3)
    if len(ifront3[0]): ifrontl[3][i] = ifront3[0][0]; ifrontu[3][i] = ifront3[0][-1]
    ifront4 = np.nonzero(inside[i]>=1e4)
    if len(ifront4[0]): ifrontl[4][i] = ifront4[0][0]; ifrontu[4][i] = ifront4[0][-1]
    ifront5 = np.nonzero(inside[i]>=1e5)
    if len(ifront5[0]): ifrontl[5][i] = ifront5[0][0]; ifrontu[5][i] = ifront5[0][-1]
    pass
for i in range(6): ax_iwf.fill_between(time, ifrontl[i]*dx, ifrontu[i]*dx, label='$n = 10^{}$ mm$^{{-2}}$'.format(i))
ax_iwf.xaxis.set_major_locator(plt.MultipleLocator(12)) # Format time axis ticks
ax_iwf.set_xlabel('Time (hrs)')
ax_iwf.set_ylabel('Distance from catheter tip (mm)')
ax_iwf.set_title('Wavefront inside')
ax_iwf.legend(loc='upper right')

plt.show()