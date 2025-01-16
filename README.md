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
