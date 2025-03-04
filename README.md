# README

[![DOI](https://zenodo.org/badge/814253783.svg)](https://doi.org/10.5281/zenodo.14674826)

This project provides a set of classes, methods and a program in order to 
simulate the spread of a bacterial infection within a urinary catheter.
The model is broken down into 3 parts.
1. A Fisher wave equation is used to model the spread of bacteria up the outside
of the catheter tubing.
2. The bladder is modelled as a well-mixed volume, with an influx of bacteria
from the top of the outer catheter.
3. The inside of the catheter is modelled as a bacteria containing fluid 
flowing down a pipe under gravity. The adhesion of bacteria to 
the surface is modelled as diffusion to an absorbing surface. On the surface 
of the catheter the behaviour of the bacteria are again modelled by a Fisher
wave equation, this time with a source term.

## Requirements

This project has been tested on:
- WSL2 Ubuntu 22.04.3
- Linux RHEL 7.8

Makefile provided for GNU gcc compiler.

### BasicCatheter dependencies:
- C++ standard

### Sensitivity analysis dependencies:
- Python 3
```
numpy
pandas
matplotlib
seaborn
SALib
```

## Installation
```
git clone https://github.com/freyabull/biophysical-catheter-model.git
cd biophysical-catheter-model
make
```

## Running
### BasicCatheter demo
```
cd BasicCatheter/BasicCatheter
./Catheter
```
The file `params.txt` contains the name of the file results will be written to (default: `results.csv`), followed by the simulation parameters.

### Sensitivity analysis
The file `sensitivity_analysis/sensitivity_analysis.ipynb` contains an interactive (Jupyter) notebook. Open the file and `Run All` to generate a sensitivity analysis of the provided dataset.

### Data analysis
The file `sensitivity_analysis/single_data.py` is an example python script for analysing the `results.csv` file generated by `BasicCatheter` to find the maximum bladder density and the time taken to attain 0.99*maximum.
