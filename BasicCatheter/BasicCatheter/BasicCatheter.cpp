// BasicCatheter.cpp : This file contains the 'main' function. Program
// execution begins and ends there.
// Simulates a urinary catheter within the bladder by:
//   1. A Fisher wave equation to model bacteria spreading up the outside of
//      the catheter.
//   2. Modelling the residual urine within the bladder as a well mixed volume.
//   3. Modelling the inside of the catheter as droplets of contaminated urine
//      rolling down a surface under gravity, with bacteria within the droplet
//      having some fixed "stickiness".
//

#include "pch.h"
#include <iostream>
#include "Catheter.h"
#include "BasicParameters.h"
#include "PDE.h"

int main()
{
    // Define variables
	double diffusivity = 1e-2; // diffusivity of bacteria on catheter, in mm^2/s
	double growth_rate1 = 1 / 3600; // Growth rate of bacteria on outside of catheter, in reproductions/s
	double carrying_capacity1 = 1e6; // Carrying capacity of bacteria on outside of catheter, per mm^2
	double growth_rate2 = 2 / 3600; // Growth rate of bacteria in bladder, in reproductions/s
	double carrying_capacity2 = 1e8; // Carrying capacity of bacteria in bladder, per mL
	double growth_rate3 = 1 / 3600; // Growth rate of bacteria on inside of catheter, in reproductions/s
	double carrying_capacity3 = 1e6; // Carrying capacity of bacteria on inside of catheter, per mm^2
	double urine_rate = 1 / 60; // Rate of generation of urine by kidneys, in mL/min
	double droplet_size = 0.05; // Size of the average droplet of urine, in mL
	double stickiness = 1; // Probability of a bacterium sticking to surface if it comes into contact
	double sump_volume = 10; // Volume of residual urine in bladder, in mL
	double catheter_length = 160; // Length of catheter, in mm
	double skin_concentration = 1e2; // Concentration of bacteria on skin (ie boundary condition)
	int x_len = 10; // Number of x steps

	Catheter myCatheter = Catheter(skin_concentration, x_len);
	BasicParameters myParam = BasicParameters(diffusivity, growth_rate1,
		carrying_capacity1, growth_rate2, carrying_capacity2, growth_rate3,
		carrying_capacity3, urine_rate, droplet_size, stickiness, sump_volume, catheter_length);

	PDE myPDE = PDE(&myParam, &myCatheter, 604800, 10, 3600);
	myPDE.solve();

}
