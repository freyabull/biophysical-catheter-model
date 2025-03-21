# README

## Raw data generated by catheter model

### Initial conditions
Raw data for four different initial/boundary conditions can be found in the folder `initial_conditions`. Each csv file contains the dynamics of bacterial colonisation originating at:
- the drainage bag (`results_bag.csv`)
- the bladder (`results_bladder.csv`)
- the skin-catheter interface (`results_skin.csv`)
- an even distribution of bacteria on the external catheter surface at insertion (`results_uniform.csv`)

### Urethral lengths
Raw data for the effect of varying urethral length can be found in the folder `urethral_lengths`. The file `parameters.txt` contains a list of parameters used for each simulation: each line $i$ in the file are the parameters for one simulation, with the corresponding results found in `results/urethral_length_i.csv`.

### Urine production rates
Raw data for the effect of varying urine production rate can be found in the folder `urine_rates`. The file `parameters.txt` contains a list of parameters used for each simulation: each line $i$ in the file are the parameters for one simulation, with the corresponding results found in `results/urine_rate_i.csv`.


## Interpretting csv files
The first two lines of the csv file contain the simulation metadata. The third line will be empty. Starting at line 4, for each output timestep $t$:
- line $4t$ will contain the bacteria per $mm^2$ on the external catheter surface (for each output position $x$)
- line $4t+1$ will contain the bacteria per $mm^3$ in the residual urine in the bladder
- line $4t+2$ will contain the bacteria per $mm^2$ on the internal catheter surface (for each output position $x$)
- line $4t+3$ will contain the mean bacteria per $mm^3$ in the outflowing urine, estimated by trapezoidal integration, and the error on this estimate.
