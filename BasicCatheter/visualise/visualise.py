""" Module to make graphs to visualise the catheter simulation. """

import matplotlib.pyplot as plt 
import pandas as pd

data = pd.read_csv('../BasicCatheter/results.csv', index_col=0, header=None).T
print(data)
