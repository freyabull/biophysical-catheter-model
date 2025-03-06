import numpy as np
import pandas as pd

i = INDEX
file_name = FILENAME
output_file = 'output_'+str(i)+'.txt'
proceed = True
try:
    df = pd.read_csv(file_name, index_col=0, header=None, skiprows=3)
    info = pd.read_csv(file_name,  nrows=1)
except:
    print('unable to open file ', file_name, ' at index ', i)
    proceed = False
if df.empty: proceed = False

if proceed:
    t_len = int(info['simulation length'].iloc[0]/info['print interval'].iloc[0]) # number of time steps
    dt = float(info['print interval'].iloc[0]/3600) # print interval (hrs)

    bladder = np.array([df.iloc[4*j+1][1] for j in range(0,t_len)])

    # Find maximum bladder infection & time to infection
    ss_bac = np.max(bladder)
    time_b = dt*np.searchsorted(bladder,0.99*ss_bac, side='right')
    final_o = df.iloc[4*(t_len-1),1]
    final_i = df.iloc[4*(t_len-1)+2,-1]
    valid = final_i>0.99*final_o # Check simulation has attained steady state

    output_array = [i, ss_bac, time_b, valid]
    np.savetxt(output_file, output_array)