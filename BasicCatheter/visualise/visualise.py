""" Module to make graphs to visualise the catheter simulation. """

import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd

plt.rc('text', usetex=True) # Enable latex

df = pd.read_csv('../BasicCatheter/results.csv', index_col=0, header=None, skiprows=3)
info = pd.read_csv('../BasicCatheter/results.csv',  nrows=1)
print(df)
print(info)

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
po_min = np.zeros(t_len)
po_max = np.zeros(t_len)
for i in range(1,t_len):
    po_range = np.nonzero(outside[i]>=(outside[i][int(peak_outside[i])]*0.99))
    po_min[i] = po_range[0][0]
    po_max[i] = po_range[0][-1]

# Find inside peaks
peak_inside = np.argmax(inside, axis=1)
pi_min = np.zeros(t_len)
pi_max = np.zeros(t_len)
for i in range(1,t_len):
    pi_range = np.nonzero(inside[i]>=(inside[i][int(peak_inside[i])]*0.99))
    pi_min[i] = pi_range[0][0]
    pi_max[i] = pi_range[0][-1]

# Find wavefronts outside
ofront1l = np.empty(t_len)
ofront1u = np.empty(t_len)
ofront2l = np.empty(t_len)
ofront2u = np.empty(t_len)
ofront3l = np.empty(t_len)
ofront3u = np.empty(t_len)
ofront4l = np.empty(t_len)
ofront4u = np.empty(t_len)
ofront5l = np.empty(t_len)
ofront5u = np.empty(t_len)
ofront0l = np.empty(t_len)
ofront0u = np.empty(t_len)
ofront1l[:]= np.nan
ofront1u[:]= np.nan
ofront2l[:]= np.nan
ofront2u[:]= np.nan
ofront3l[:]= np.nan
ofront3u[:]= np.nan
ofront4l[:]= np.nan
ofront4u[:]= np.nan
ofront5l[:]= np.nan
ofront5u[:]= np.nan
ofront0l[:]= np.nan
ofront0u[:]= np.nan
for i in range(1,t_len):
    ofront1 = np.nonzero(outside[i]>10e1)
    if len(ofront1[0]):
        ofront1l[i] = ofront1[0][0]
        ofront1u[i] = ofront1[0][-1]
    ofront2 = np.nonzero(outside[i]>10e2)
    if len(ofront2[0]):
        ofront2l[i] = ofront2[0][0]
        ofront2u[i] = ofront2[0][-1]
    ofront3 = np.nonzero(outside[i]>10e3)
    if len(ofront3[0]):
        ofront3l[i] = ofront3[0][0]
        ofront3u[i] = ofront3[0][-1]
    ofront4 = np.nonzero(outside[i]>10e4)
    if len(ofront4[0]):
        ofront4l[i] = ofront4[0][0]
        ofront4u[i] = ofront4[0][-1]
    ofront5 = np.nonzero(outside[i]>10e5)
    if len(ofront5[0]):
        ofront5l[i] = ofront5[0][0]
        ofront5u[i] = ofront5[0][-1]
    ofront0 = np.nonzero(outside[i]>10e0)
    if len(ofront0[0]):
        ofront0l[i] = ofront0[0][0]
        ofront0u[i] = ofront0[0][-1]

# Find wavefronts inside
ifront1l = np.empty(t_len)
ifront1u = np.empty(t_len)
ifront2l = np.empty(t_len)
ifront2u = np.empty(t_len)
ifront3l = np.empty(t_len)
ifront3u = np.empty(t_len)
ifront4l = np.empty(t_len)
ifront4u = np.empty(t_len)
ifront5l = np.empty(t_len)
ifront5u = np.empty(t_len)
ifront0l = np.empty(t_len)
ifront0u = np.empty(t_len)
ifront1l[:]= np.nan
ifront1u[:]= np.nan
ifront2l[:]= np.nan
ifront2u[:]= np.nan
ifront3l[:]= np.nan
ifront3u[:]= np.nan
ifront4l[:]= np.nan
ifront4u[:]= np.nan
ifront5l[:]= np.nan
ifront5u[:]= np.nan
ifront0l[:]= np.nan
ifront0u[:]= np.nan
for i in range(1,t_len):
    ifront1 = np.nonzero(inside[i]>10e1)
    if len(ifront1[0]):
        ifront1l[i] = ifront1[0][0]
        ifront1u[i] = ifront1[0][-1]
    ifront2 = np.nonzero(inside[i]>10e2)
    if len(ifront2[0]):
        ifront2l[i] = ifront2[0][0]
        ifront2u[i] = ifront2[0][-1]
    ifront3 = np.nonzero(inside[i]>10e3)
    if len(ifront3[0]):
        ifront3l[i] = ifront3[0][0]
        ifront3u[i] = ifront3[0][-1]
    ifront4 = np.nonzero(inside[i]>10e4)
    if len(ifront4[0]):
        ifront4l[i] = ifront4[0][0]
        ifront4u[i] = ifront4[0][-1]
    ifront5 = np.nonzero(inside[i]>10e5)
    if len(ifront5[0]):
        ifront5l[i] = ifront5[0][0]
        ifront5u[i] = ifront5[0][-1]
    ifront0 = np.nonzero(inside[i]>10e0)
    if len(ifront0[0]):
        ifront0l[i] = ifront0[0][0]
        ifront0u[i] = ifront0[0][-1]


fig, ((axb, axo, axi, ax), (axow, axiw, axwo, axwi)) = plt.subplots(2,4)

# Plot bladder concentration against time
axb.plot(time[6:], bladder[6:])
axb.xaxis.set_major_locator(plt.MultipleLocator(12))
axb.set_yscale('log')
#axb.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
axb.set_xlabel('Time (hrs)')
axb.set_ylabel('Bacterial concentration (mm$^{-3}$)')
axb.set_title('Bladder')

# Plot wave profile outside
axo.plot(x,outside[16])
axo.plot(x,outside[18])
axo.plot(x,outside[20])
axo.plot(x,outside[24])
axo.plot(x,outside[47])
axo.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
axo.set_xlabel('Distance from catheter base (mm)')
axo.set_ylabel('Bacterial concentration (mm$^{-2}$)')
axo.set_title('Outside')

# Plot wave profile inside
axi.plot(x,inside[22])
axi.plot(x,inside[24])
axi.plot(x,inside[26])
axi.plot(x,inside[28])
axi.plot(x,inside[47])
axi.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
axi.set_xlabel('Distance from catheter tip (mm)')
axi.set_ylabel('Bacterial concentration (mm$^{-2}$)')
axi.set_title('Inside large time')

#Plot wavepeak outside
axow.plot(time[1:], peak_outside[1:]*dx)
axow.fill_between(time, po_min*dx, po_max*dx, alpha=0.5)
axow.xaxis.set_major_locator(plt.MultipleLocator(12))
axow.set_xlabel('Time (hrs)')
axow.set_ylabel('Distance from catheter base (mm)')
axow.set_title('Position of extraluminal peak concentration')

#Plot wavepeak inside
axiw.plot(time[1:], peak_inside[1:]*dx)
axiw.fill_between(time, pi_min*dx, pi_max*dx, alpha=0.5)
axiw.xaxis.set_major_locator(plt.MultipleLocator(12))
axiw.set_xlabel('Time (hrs)')
axiw.set_ylabel('Distance from catheter tip (mm)')
axiw.set_title('Position of intraluminal peak concentration')

#Plot inside small time
ax.plot(x,inside[7])
ax.plot(x,inside[9])
ax.plot(x,inside[11])
ax.plot(x,inside[13])
ax.plot(x,inside[15])
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax.set_xlabel('Distance from catheter tip (mm)')
ax.set_ylabel('Bacterial concentration (mm$^{-2}$)')
ax.set_title('Inside small time')

#Plot wavefront outside
axwo.fill_between(time, ofront1l*dx, ofront1u*dx, alpha=0.2, color='b', label='10e1')
axwo.fill_between(time, ofront2l*dx, ofront2u*dx, alpha=0.2, color='b', label='10e2')
axwo.fill_between(time, ofront3l*dx, ofront3u*dx, alpha=0.2, color='b', label='10e3')
axwo.fill_between(time, ofront4l*dx, ofront4u*dx, alpha=0.2, color='b', label='10e4')
axwo.fill_between(time, ofront5l*dx, ofront5u*dx, alpha=0.2, color='b', label='10e5')
axwo.fill_between(time, ofront0l*dx, ofront0u*dx, alpha=0.2, color='b', label='10e0')
axwo.xaxis.set_major_locator(plt.MultipleLocator(12))
axwo.set_xlabel('Time (hrs)')
axwo.set_ylabel('Distance from catheter base (mm)')
axwo.set_title('Wavefront outside')

#Plot wavefront inside
axwi.fill_between(time, ifront1l*dx, ifront1u*dx, alpha=0.2, color='b', label='10e1')
axwi.fill_between(time, ifront2l*dx, ifront2u*dx, alpha=0.2, color='b', label='10e2')
axwi.fill_between(time, ifront3l*dx, ifront3u*dx, alpha=0.2, color='b', label='10e3')
axwi.fill_between(time, ifront4l*dx, ifront4u*dx, alpha=0.2, color='b', label='10e4')
axwi.fill_between(time, ifront5l*dx, ifront5u*dx, alpha=0.2, color='b', label='10e5')
axwi.fill_between(time, ifront0l*dx, ifront0u*dx, alpha=0.2, color='b', label='10e0')
axwi.xaxis.set_major_locator(plt.MultipleLocator(12))
axwi.set_xlabel('Time (hrs)')
axwi.set_ylabel('Distance from catheter tip (mm)')
axwi.set_title('Wavefront inside')

plt.show()