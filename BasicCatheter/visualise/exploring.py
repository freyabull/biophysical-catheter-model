""" Module to make graphs from multiple catheter simulations. """

import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import seaborn as sns
import os

# Location of results files
folder = '../ExploringCatheter/results/urine_rates'
files = os.listdir(folder)

# Arrays for results
N = len(files)
urine_rates = np.zeros(N)
times = np.zeros(N)

infection_threshold = 1e5 # Clinical definition of infection occuring

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
    outflow = np.zeros(t_len)
    for j in range(0, t_len):
        outflow[j] = df.iloc[4*j+3][1]

    index = np.searchsorted(outflow,infection_threshold, side='right')
    if index==len(outflow):
        times[i] = np.nan
    else :
        times[i] = dt*index
 
plt.plot(urine_rates, times) # Graph timescale to clinical infection at outlet against urine rate
plt.xlabel("Urine flow rate mm/s")
plt.ylabel("Time to infection s")
plt.show()