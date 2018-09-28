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
#include <stdlib.h>
#include <fstream>
#include "Catheter.h"
#include "BasicParameters.h"
#include "PDE.h"
#include "../GetOptLibrary/getopt.h"

void showHelp() {
	std::cout << "This is the help for BasicCatheter." << std::endl;
	std::cout << "\t -d \t Diffusivity of bacteria (mm^2/s) \t\t Default = 1E-2" << std::endl;
	std::cout << "\t -u \t Rate of generation of urine (per min) \t\t Default = 1" << std::endl;
	std::cout << "\t -s \t Droplet size (mL) \t\t\t\t Default = 0.05" << std::endl;
	std::cout << "\t -p \t Probability of bacteria sticking \t\t Default = 1" << std::endl;
	std::cout << "\t -v \t Volume of residual urine (mL) \t\t\t Default = 10" << std::endl;
	std::cout << "\t -l \t Catheter length (mm) \t\t\t\t Default = 160" << std::endl;
	std::cout << "\t -b \t Skin concentration of bacteria (per mm^2) \t Default = 1E2" << std::endl;
	std::cout << "\t -c \t Type of initial condition \t\t\t Default = ? \t Choose from:" << std::endl;
	std::cout << "\t -f \t Output file \t\t\t\t\t Default = results.csv" << std::endl;
}

int main(int argc, char **argv)
{
    // Define variables
	double diffusivity = 1e-2; // diffusivity of bacteria on catheter, in mm^2/s
	double growth_rate1 = 1.0 / 3600.0; // Growth rate of bacteria on outside of catheter, in reproductions/s
	double carrying_capacity1 = 1e6; // Carrying capacity of bacteria on outside of catheter, per mm^2
	double growth_rate2 = 2.0 / 3600.0; // Growth rate of bacteria in bladder, in reproductions/s
	double carrying_capacity2 = 1e8; // Carrying capacity of bacteria in bladder, per mL
	double growth_rate3 = 1.0 / 3600.0; // Growth rate of bacteria on inside of catheter, in reproductions/s
	double carrying_capacity3 = 1e6; // Carrying capacity of bacteria on inside of catheter, per mm^2
	double urine_rate = 1.0 / 60.0; // Rate of generation of urine by kidneys, in mL/min
	double droplet_size = 0.05; // Size of the average droplet of urine, in mL
	double stickiness = 1.0; // Probability of a bacterium sticking to surface if it comes into contact
	double sump_volume = 10.0; // Volume of residual urine in bladder, in mL
	double catheter_length = 160.0; // Length of catheter, in mm
	double skin_concentration = 1e2; // Concentration of bacteria on skin (ie boundary condition)
	char *initial_condition; // Type of initial conditions
	int x_len = 10; // Number of x steps
	const char *file_name = "results.csv"; // File for output to be written to


	int c;
	opterr = 0;
	char *endptr;

	while((c = getopt(argc, argv, "hd:u:s:p:v:l:b:c:f:")) != -1)
		switch (c)
		{
		case 'h':
			showHelp();
			break;
		case 'd':
			diffusivity = strtod(optarg, &endptr);
			break;
		case 'u':
			urine_rate = strtod(optarg, &endptr);
			break;
		case 's':
			droplet_size = strtod(optarg, &endptr);
			break;
		case 'p':
			stickiness = strtod(optarg, &endptr);
			break;
		case 'v':
			sump_volume = strtod(optarg, &endptr);
			break;
		case 'l':
			catheter_length = strtod(optarg, &endptr);
			break;
		case 'b':
			skin_concentration = strtod(optarg, &endptr);
			break;
		case 'c':
			initial_condition = optarg;
			break;
		case 'f':
			file_name = optarg;
			break;
		case '?':
			std::cout << "Invalid option " << optopt << std::endl;
		default:
			abort();
		}

	Catheter myCatheter = Catheter(skin_concentration, x_len);
	BasicParameters myParam = BasicParameters(diffusivity, growth_rate1,
		carrying_capacity1, growth_rate2, carrying_capacity2, growth_rate3,
		carrying_capacity3, urine_rate, droplet_size, stickiness, sump_volume, catheter_length);

	PDE myPDE = PDE(&myParam, &myCatheter, 604800, 1, 3600);
	std::ofstream results_file; // File for output to be written to 
	results_file.open(file_name, std::ios::trunc);
	if (!results_file.is_open())
	{
		std::cout << "Unable to open file " << file_name << std::endl;
		std::cout << "Writing direct to console instead" << std::endl;
		myPDE.solve();
	}
	else
	{ 
		std::cout << "Successfully opened " << file_name << std::endl;
		myPDE.solve(results_file);
	}

	results_file.close();

}
