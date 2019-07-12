""" Module to make graphs to visualise the catheter simulation. """

import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import seaborn as sns
from textwrap import wrap

plt.rc('text', usetex=True) # Enable latex

df = pd.read_csv('../BasicCatheter/results.csv', index_col=0, header=None, skiprows=3)
info = pd.read_csv('../BasicCatheter/results.csv',  nrows=1)
print(df)
print(info)

# Find reasonable minimum bladder concentration
min_bladder = float(1.0/(1000.0*info['sump volume']))

# Extract time series
t_len = int(info['simulation length']/info['print step']) # number of time steps
dt = float(info['print step']/3600) # print interval
time = dt*np.array(range(0, t_len))

# Extract x series
x_len = int(info['num of x steps'])
dx = float(info['catheter length']/(x_len-1))
x = dx*np.array(range(0,x_len))

# Extract bladder time series
bladder = np.zeros(t_len)
for i in range(0, t_len):
    bladder[i] = df.iloc[3*i+1][1]

# Extract outside
outside = np.zeros((t_len,df.shape[1]))
for i in range(0,t_len):
    outside[i] = df.iloc[3*i]

# Extract inside
inside = np.zeros((t_len,df.shape[1]))
for i in range(0,t_len):
    inside[i] = df.iloc[3*i+2]

# Find outside peaks
peak_outside = np.argmax(outside, axis=1)
opeaks = np.array([outside[i][peak_outside[i]] for i in range(1,t_len)])
otimes = np.unique([np.argmax(opeaks>opeaks[-1]*0.2*i) for i in range(1,5)])
po_min = np.empty(t_len)
po_max = np.empty(t_len)
po_min[:] = np.nan
po_max[:] = np.nan
for i in range(1,t_len):
    po_range = np.nonzero(outside[i]>=(outside[i][int(peak_outside[i])]*0.99))
    if len(po_range[0]):
        po_min[i] = po_range[0][0]
        po_max[i] = po_range[0][-1]

# Find inside peaks
peak_inside = np.argmax(inside, axis=1)
ipeaks = np.array([inside[i][peak_inside[i]] for i in range(1,t_len)])
itimes = np.unique([np.argmax(ipeaks>ipeaks[-1]*0.2*i) for i in range(1,5)])
istimes = np.unique([np.argmax(ipeaks>(20*i)) for i in range(1,5)])
pi_min = np.empty(t_len)
pi_max = np.empty(t_len)
pi_min[:] = np.nan
pi_max[:] = np.nan
for i in range(1,t_len):
    pi_range = np.nonzero(inside[i]>=(inside[i][int(peak_inside[i])]*0.99))
    if len(pi_range[0]):
        pi_min[i] = pi_range[0][0]
        pi_max[i] = pi_range[0][-1]

# Find wavefronts outside
ofrontl = np.empty([6,t_len])
ofrontu = np.empty([6,t_len])
ofrontl[:][:] = np.nan
ofrontu[:][:] = np.nan
for i in range(1,t_len):
    ofront0 = np.nonzero(outside[i]>=10e0)
    if len(ofront0[0]):
        ofrontl[0][i] = ofront0[0][0]
        ofrontu[0][i] = ofront0[0][-1]
    ofront1 = np.nonzero(outside[i]>=10e1)
    if len(ofront1[0]):
        ofrontl[1][i] = ofront1[0][0]
        ofrontu[1][i] = ofront1[0][-1]
    ofront2 = np.nonzero(outside[i]>=10e2)
    if len(ofront2[0]):
        ofrontl[2][i] = ofront2[0][0]
        ofrontu[2][i] = ofront2[0][-1]
    ofront3 = np.nonzero(outside[i]>=10e3)
    if len(ofront3[0]):
        ofrontl[3][i] = ofront3[0][0]
        ofrontu[3][i] = ofront3[0][-1]
    ofront4 = np.nonzero(outside[i]>=10e4)
    if len(ofront4[0]):
        ofrontl[4][i] = ofront4[0][0]
        ofrontu[4][i] = ofront4[0][-1]
    ofront5 = np.nonzero(outside[i]>=10e5)
    if len(ofront5[0]):
        ofrontl[5][i] = ofront5[0][0]
        ofrontu[5][i] = ofront5[0][-1]

# Find wavefronts inside
ifrontl = np.empty([6,t_len])
ifrontu = np.empty([6,t_len])
ifrontl[:][:] = np.nan
ifrontu[:][:] = np.nan
for i in range(1,t_len):
    ifront0 = np.nonzero(inside[i]>=10e0)
    if len(ifront0[0]):
        ifrontl[0][i] = ifront0[0][0]
        ifrontu[0][i] = ifront0[0][-1]
    ifront1 = np.nonzero(inside[i]>=10e1)
    if len(ifront1[0]):
        ifrontl[1][i] = ifront1[0][0]
        ifrontu[1][i] = ifront1[0][-1]
    ifront2 = np.nonzero(inside[i]>=10e2)
    if len(ifront2[0]):
        ifrontl[2][i] = ifront2[0][0]
        ifrontu[2][i] = ifront2[0][-1]
    ifront3 = np.nonzero(inside[i]>=10e3)
    if len(ifront3[0]):
        ifrontl[3][i] = ifront3[0][0]
        ifrontu[3][i] = ifront3[0][-1]
    ifront4 = np.nonzero(inside[i]>=10e4)
    if len(ifront4[0]):
        ifrontl[4][i] = ifront4[0][0]
        ifrontu[4][i] = ifront4[0][-1]
    ifront5 = np.nonzero(inside[i]>=10e5)
    if len(ifront5[0]):
        ifrontl[5][i] = ifront5[0][0]
        ifrontu[5][i] = ifront5[0][-1]

sns.set_palette(sns.dark_palette("sky blue", n_colors=8, input="xkcd", reverse=True))
fig, ((axb, axo, axi, ax), (axow, axiw, axwo, axwi)) = plt.subplots(2,4)
fig.suptitle('\n'.join(wrap('Model parameters: {}'.format([(i,info[i][0]) for i in info]),300)), fontsize=10)

# Plot bladder concentration against time
bladder = bladder*(bladder>min_bladder) # Remove unphysical small values
axb.plot(time, bladder)
axb.xaxis.set_major_locator(plt.MultipleLocator(12))
axb.set_yscale('log')
axb.set_xlabel('Time (hrs)')
axb.set_ylabel('Bacterial concentration (mm$^{-3}$)')
axb.set_title('Bladder')

# Plot wave profile outside
for i in otimes:
    axo.plot(x,outside[i], label='$t = {}$ hrs'.format(time[i]))
axo.plot(x,outside[-1], label='$t = {}$ hrs'.format(time[-1]))
axo.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
axo.set_xlabel('Distance from catheter base (mm)')
axo.set_ylabel('Bacterial concentration (mm$^{-2}$)')
axo.set_title('Outside')
axo.legend(loc='upper right')

# Plot wave profile inside
for i in itimes:
    axi.plot(x,inside[i], label='$t = {}$ hrs'.format(time[i]))
axi.plot(x,inside[-1],label='$t = {}$ hrs'.format(time[-1]))
axi.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
axi.set_xlabel('Distance from catheter tip (mm)')
axi.set_ylabel('Bacterial concentration (mm$^{-2}$)')
axi.set_title('Inside late time')
axi.legend(loc='upper right')

#Plot wavepeak outside
axow.plot(time[1:], peak_outside[1:]*dx)
axow.fill_between(time, po_min*dx, po_max*dx, alpha=0.2)
axow.xaxis.set_major_locator(plt.MultipleLocator(12))
axow.set_xlabel('Time (hrs)')
axow.set_ylabel('Distance from catheter base (mm)')
axow.set_title('Position of extraluminal peak concentration')

#Plot wavepeak inside
axiw.plot(time[1:], peak_inside[1:]*dx)
axiw.fill_between(time, pi_min*dx, pi_max*dx, alpha=0.2)
axiw.xaxis.set_major_locator(plt.MultipleLocator(12))
axiw.set_xlabel('Time (hrs)')
axiw.set_ylabel('Distance from catheter tip (mm)')
axiw.set_title('Position of intraluminal peak concentration')

#Plot inside small time
for i in istimes:
    ax.plot(x,inside[i], label='$t = {}$ hrs'.format(time[i]))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax.set_xlabel('Distance from catheter tip (mm)')
ax.set_ylabel('Bacterial concentration (mm$^{-2}$)')
ax.set_title('Inside early time')
ax.legend(loc='upper right')

#Plot wavefront outside
for i in range(6):
    axwo.fill_between(time, ofrontl[i]*dx, ofrontu[i]*dx, label='$n = 10^{}$ mm$^{{-2}}$'.format(i))
axwo.xaxis.set_major_locator(plt.MultipleLocator(12))
axwo.set_xlabel('Time (hrs)')
axwo.set_ylabel('Distance from catheter base (mm)')
axwo.set_title('Wavefront outside')
axwo.legend(loc='upper right')



#Plot wavefront inside
for i in range(6):
    axwi.fill_between(time, ifrontl[i]*dx, ifrontu[i]*dx, label='$n = 10^{}$ mm$^{{-2}}$'.format(i))
axwi.xaxis.set_major_locator(plt.MultipleLocator(12))
axwi.set_xlabel('Time (hrs)')
axwi.set_ylabel('Distance from catheter tip (mm)')
axwi.set_title('Wavefront inside')
axwi.legend(loc='upper right')


plt.show()