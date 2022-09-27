""" Module to make graphs from multiple catheter simulations. """

import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import seaborn as sns
import os
from matplotlib.font_manager import FontProperties
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker

# Graphics settings
plt.rcParams["figure.figsize"] = (3.1,3.1)
plt.rcParams['xtick.major.pad']='4'
plt.rcParams['ytick.major.pad']='4'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = '8'
font = FontProperties()
font.set_name('serif')
font.set_size(8)
colours = np.array([(68,120,33), (141,211,95), (255,179,128), (255,221,85), (179,179,179)])/255
palette = sns.color_palette(colours)
cmappy = LinearSegmentedColormap.from_list("whitepink",[palette[2],'white'], N=256)
cmappy2 = LinearSegmentedColormap.from_list("greenwhite",['white',palette[1]], N=256)
sns.set_palette(palette)
sns.set_style("ticks")


# Process data
reimport_data = False
if reimport_data:
    # Location of results files
    folder_prefix = '../ExploringCatheter/results'
    folders = ['d4', 'd35', 'd5', 'd36', 'd6', 'd37', 'd7', 'd38', 'd8']
    files = [os.path.join(path, file) for folder in folders for path, subdirs, files in os.walk(os.path.join(folder_prefix, folder)) for file in files]

    #  Threshold value for blockage
    blockage_thresh = 1e7

    # Arrays for results
    N = len(files)
    urine_rates = np.zeros(N)
    diffusion_coeffs = np.zeros(N)
    max_infect = np.zeros(N)
    time_to_block = np.zeros(N)

    # Open files and read in data
    for i in range(N):
        filename = files[i]
        print(filename)
        # Read data
        df = pd.read_csv(filename, index_col=0, header=None, skiprows=3)
        info = pd.read_csv(filename,  nrows=1)
        if df.empty: break

        urine_rates[i] = info['urine rate']
        diffusion_coeffs[i] = info['surface diffusivity']
        t_len = int(info['simulation length']/info['print interval']) # number of time steps
        dt = float(info['print interval']/3600) # print interval (hrs)

        bladder = np.zeros(t_len)
        inside = np.zeros((t_len,df.shape[1]))
        outflow = np.zeros(t_len)
        maxy = np.zeros(t_len)
        for j in range(0, t_len):
            bladder[j] = df.iloc[4*j+1][1]
            inside[j] = df.iloc[4*j+2]
            maxy[j] = max(inside[j])

        # Find maximum bladder infection & time to infection
        max_infect[i] = np.max(bladder)

        # Find times at which inside blocks
        index = np.searchsorted(maxy,blockage_thresh, side='right')
        time_to_block[i] = dt*index/24 if index!=t_len else np.nan

    # Store results in a dataframe, save as a csv for later use
    results_frame = pd.DataFrame({'urine rate' : urine_rates, 
        'diffusion coefficients' : diffusion_coeffs, 'bladder steady state' : max_infect,
        'time to blockage' : time_to_block})
    results_frame.to_csv('parameter_space.csv',index=False)   
else:
    # Load dataframe of results
    results_frame = pd.read_csv('parameter_space.csv')


# Make graphs
# Set up data frames
results_frame = results_frame.set_index(['diffusion coefficients', 'urine rate'])
heatmap_frame = results_frame['time to blockage'].reset_index().pivot(
    columns='diffusion coefficients',index='urine rate',values='time to blockage')
contourmap_frame = results_frame['bladder steady state'].reset_index().pivot(
    columns='diffusion coefficients',index='urine rate',values='bladder steady state')
print(results_frame)

fig, ax = plt.subplots()

# Create heat map for time to blockage
Xh=heatmap_frame.columns.values
Yh=heatmap_frame.index.values
Zh=heatmap_frame.values
hm = ax.contourf(Xh,Yh,Zh,levels=6,cmap=cmappy)
# Create a colour bar
cbar = fig.colorbar(hm)
cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontproperties=font)
locator, formatter = cbar._get_ticker_locator_formatter()
locator.set_params(nbins=4)
cbar.update_ticks()
plt.text(0.006,11.5,'Time to blockage (days)',rotation=90,fontproperties=font)
ax.set_xscale('log')

# Create contour map for bacterial density
X=contourmap_frame.columns.values
Y=contourmap_frame.index.values
Z=contourmap_frame.values
cm=ax.contour(X,Y,Z,locator=plt.LogLocator(),colors='k')
# Add contour labels
fmt = ticker.LogFormatterMathtext()
fmt.create_dummy_axis()
manual_locations = [] # Need to assign suitable location values
paths = np.sum([True for line in cm.collections if line.get_paths()])
begin = contourmap_frame.index[0] # smallest value on y axis
step = (contourmap_frame.index[-1]-begin)/(paths+2) # step between placements for clabel
for line in cm.collections:
    for path in line.get_paths():
        v = path.vertices
        y = begin+(len(manual_locations)+1)*step
        x_index = np.argmin(np.abs(v[:,1]-y))
        manual_locations.append((v[x_index][0],y))
labels = ax.clabel(cm, inline=True, inline_spacing=3, fmt=fmt, 
                   manual=manual_locations) # Add contour labels
for l in labels:
    l.set_rotation(0)
    l.set_font_properties(font)

# General figure properties
plt.yticks(fontproperties=font)
plt.xticks(fontproperties=font)
plt.xlabel('Diffusion coefficient (mm)',fontproperties=font, labelpad=1.5)
plt.ylabel('Urine flow rate (mm$^3$s$^{-1}$)',fontproperties=font, labelpad=2)
plt.locator_params(axis='y',tight=True, nbins=4)
plt.locator_params(axis='x',tight=True, numticks=4)
plt.tight_layout(rect=[-0.04,-0.045,1.046,1.035])
plt.savefig('ps_diagram.pdf')

#Separate plots
fig2, ax2 = plt.subplots()
# Create heat map for time to blockage
hm2 = ax2.contourf(Xh,Yh,Zh,levels=6,cmap=cmappy)
cbar2 = fig2.colorbar(hm2)
cbar2.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontproperties=font)
locator, formatter = cbar2._get_ticker_locator_formatter()
locator.set_params(nbins=4)
cbar2.update_ticks()
ax2.set_xscale('log')
plt.text(0.004,11.5,'Time to thick biofilm (days)',rotation=90,fontproperties=font)
plt.yticks(fontproperties=font)
plt.xticks(fontproperties=font)
plt.xlabel('Diffusion coefficient (mm)',fontproperties=font, labelpad=1.5)
plt.ylabel('Urine flow rate (mm$^3$s$^{-1}$)',fontproperties=font, labelpad=2)
plt.locator_params(axis='y',tight=True, nbins=4)
plt.locator_params(axis='x',tight=True, numticks=4)
plt.tight_layout(rect=[-0.04,-0.045,1.046,1.035])
plt.savefig('ps_diagram_t2b.pdf')
plt.savefig('ps_diagram_t2b.svg')

#Separate plots
fig3, ax3 = plt.subplots()
# Create contour map for bacterial density
hm3=ax3.contourf(X,Y,Z,locator=plt.LogLocator(),cmap=cmappy2)
cbar3 = fig3.colorbar(hm3)
cbar3.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontproperties=font)
locator, formatter = cbar3._get_ticker_locator_formatter()
locator.set_params(nbins=4)
cbar3.update_ticks()
ax3.set_xscale('log')
plt.text(0.0017,11.7,'Max bladder density (mm$^{-3}$)',rotation=90,fontproperties=font)
plt.yticks(fontproperties=font)
plt.xticks(fontproperties=font)
plt.xlabel('Diffusion coefficient (mm)',fontproperties=font, labelpad=1.5)
plt.ylabel('Urine flow rate (mm$^3$s$^{-1}$)',fontproperties=font, labelpad=2)
plt.locator_params(axis='y',tight=True, nbins=4)
plt.locator_params(axis='x',tight=True, numticks=4)
plt.tight_layout(rect=[-0.04,-0.045,1.048,1.035])
plt.savefig('ps_diagram_md.pdf')
plt.savefig('ps_diagram_md.svg')

plt.show()
