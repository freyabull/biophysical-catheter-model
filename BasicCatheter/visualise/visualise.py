""" Module to make graphs to visualise the catheter simulation. """

import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd

plt.rc('text', usetex=True) # Enable latex

df = pd.read_csv('../BasicCatheter/results.csv', index_col=0, header=None, skiprows=3)
info = pd.read_csv('../BasicCatheter/results.csv',  nrows=1)
print(df)
print(info)

# Extract bladder time series
time = np.array(range(0, 24))
bladder = np.zeros(24)
for i in range(0, 24):
    bladder[i] = df.iloc[3*i+1][1]

# Extract outside
outside = np.zeros((24,df.shape[1]))
for i in range(0,24):
    outside[i] = df.iloc[3*i]

# Extract inside
inside = np.zeros((24,df.shape[1]))
for i in range(0,24):
    inside[i] = df.iloc[3*i+2]

# Find outside peaks
peak_outside = np.argmax(outside, axis=1)

# Find inside peaks
peak_inside = np.argmax(inside, axis=1)

fig, ((axb, axo, axi), (axow, axiw, ax)) = plt.subplots(2,3)

# Plot bladder concentration against time
axb.plot(time, bladder)
axb.set_xlabel('Time (hrs)')
axb.set_ylabel('Bacterial concentration (mm$^{-3}$)')
axb.set_title('Bladder')

# Plot wave profile outside
axo.plot(outside[14])
axo.plot(outside[17])
axo.plot(outside[20])
axo.plot(outside[23])
axo.set_xlabel('\% Distance up catheter')
axo.set_ylabel('Bacterial concentration (mm$^{-2}$)')
axo.set_title('Outside')

# Plot wave profile inside
axi.plot(inside[20])
axi.plot(inside[21])
axi.plot(inside[22])
axi.plot(inside[23])
axi.set_xlabel('\% Distance down catheter')
axi.set_ylabel('Bacterial concentration (mm$^{-2}$)')
axi.set_title('Inside')

#Plot wavefront outside
axow.plot(time, peak_outside)
axow.set_xlabel('Time (hrs)')
axow.set_ylabel('\% Distance up catheter')
axow.set_title('Location of peak concentration, outside')

#Plot wavefront inside
axiw.plot(time, peak_inside)
axiw.set_xlabel('Time (hrs)')
axiw.set_ylabel('\% Distance down catheter')
axiw.set_title('Location of peak concentration, inside')

#plt.show()