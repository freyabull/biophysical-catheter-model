// BasicCatheter.cpp : This file contains the 'main' function. Program
// execution begins and ends there.
// Simulates a urinary catheter within the bladder by:
//   1. A Fisher wave equation to model bacteria spreading up the outside of
//      the catheter.
//   2. Modelling the residual urine within the bladder as a well mixed volume.
//   3. Modelling the inside of the catheter as contaminated urine
//      flowing through the tube, where bacteria within a certain distance of
//      the tube walls have some probability of sticking
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
	std::cout << "\t -a \t Catheter radius (mm) \t\t\t\t Default = 1.0" << std::endl;
	std::cout << "\t -p \t Probability of bacteria sticking \t\t Default = 1" << std::endl;
	std::cout << "\t -v \t Volume of residual urine (mL) \t\t\t Default = 10" << std::endl;
	std::cout << "\t -l \t Catheter length (mm) \t\t\t\t Default = 160" << std::endl;
	std::cout << "\t -s \t Skin concentration of bacteria (per mm^2) \t Default = 1E2" << std::endl;
	std::cout << "\t -b \t Drainage bag bacterial concentration (per mL) \t Default = 1E2" << std::endl;
	std::cout << "\t -c \t Type of initial condition \t\t\t Default = 1 "
		"\n \t\t Choose from: \t 1 \t Clean catheter "
		"\n \t\t\t\t 2 \t Uniform distribution "
		"\n \t\t\t\t 3 \t Exponential distribution " 
		"\n \t\t\t\t 4 \t Infected bladder " << std::endl;
	std::cout << "\t -f \t Output file \t\t\t\t\t Default = results.csv" << std::endl;
}

int main(int argc, char **argv)
{
    // Define variables
	double diffusivity = 1e-2; // diffusivity of bacteria on catheter, in mm^2/s
	double growth_rate1 = std::log(2.0) / 3600.0; // Growth rate of bacteria on outside of catheter, per s
	double carrying_capacity1 = 1e6; // Carrying capacity of bacteria on outside of catheter, per mm^2
	double growth_rate2 = std::log(2.0) / 1800.0; // Growth rate of bacteria in bladder, per s
	double carrying_capacity2 = 1e5; // Carrying capacity of bacteria in bladder, per mm^3
	double growth_rate3 = std::log(2.0) / 3600.0; // Growth rate of bacteria on inside of catheter, per s
	double carrying_capacity3 = 1e6; // Carrying capacity of bacteria on inside of catheter, per mm^2
	double urine_rate = 1.0e3 / 60.0; // Rate of generation of urine by kidneys, in mm^3/s NOT IN USE
	double catheter_radius = 1.0; // Radius of catheter in mm
	double stickiness = 1.0; // Probability of a bacterium sticking to surface if it comes into contact
	double sump_volume = 1.0e4; // Volume of residual urine in bladder, in mm^3
	double catheter_length = 160.0; // Length of catheter, in mm
	double skin_concentration = 1e2; // Concentration of bacteria on skin (ie boundary condition)
	double bag_concentration = 0.1; // Concentration of bacteria in drainage bag (ie boundary condition)
	int initial_condition = 1; // Type of initial conditions
	int x_len = 101; // Number of x points
	const char *file_name = "results.csv"; // File for output to be written to
	double initial_skin_conc = 1e2; // Concentration of bacteria on the skin at time of catheter insertion
	double lambda = 1.0/catheter_length; // Scale factor for exponential distributions
	int simulation_length = 86400*2; // Timeframe of simulation, in s
	double dt = 1; // Time step
	int print_step = 3600; // Number of steps at which to output data

	// Get command line arguments
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
			urine_rate = strtod(optarg, &endptr)*100.0/6.0;
			break;
		case 'a':
			catheter_radius = strtod(optarg, &endptr);
			break;
		case 'p':
			stickiness = strtod(optarg, &endptr);
			break;
		case 'v':
			sump_volume = strtod(optarg, &endptr)*1000.0;
			break;
		case 'l':
			catheter_length = strtod(optarg, &endptr);
			break;
		case 's':
			skin_concentration = strtod(optarg, &endptr);
			break;
		case 'b':
			bag_concentration = strtod(optarg, &endptr)*0.001;
			break;
		case 'c':
			initial_condition = atoi(optarg);
			break;
		case 'f':
			file_name = optarg;
			break;
		case '?':
			std::cout << "Invalid option " << optopt << std::endl;
		default:
			abort();
		}

	// Set initial conditions
	Catheter myCatheter; // Current state of catheter
	std::vector<double> outside; // Concentration profile of outside of catheter
	double bladder; // Concentration within bladder
	std::vector<double> inside; // Concentration profile of inside of catheter

	switch (initial_condition) 
	{
	case 1:
		// Clean catheter
		myCatheter = Catheter(skin_concentration, bag_concentration, x_len); 
		break;
	case 2:
		// Uniform distribution
		outside = std::vector<double>(x_len, initial_skin_conc/(x_len-1));
		bladder = 0.0;
		inside = std::vector<double>(x_len, 0.0);
		myCatheter = Catheter(outside, bladder, inside, skin_concentration, bag_concentration);
		break;
	case 3:
		// Exponential distribution
		for (int i = 0; i < x_len; ++i) 
		{
			double x1 = (x_len-i-1) * catheter_length / (x_len-1);
			double x2 = (x_len-i) * catheter_length / (x_len-1);
			double value = initial_skin_conc * (exp(-lambda*x1) - exp(-lambda*x2));
			outside.push_back(value);
		}
		bladder = 0.0;
		inside = std::vector<double>(x_len, 0.0);
		myCatheter = Catheter(outside, bladder, inside, skin_concentration, bag_concentration);
		break;
	case 4:
		// Infected bladder
		outside = std::vector<double>(x_len, 0.0);
		bladder = 1e2;
		inside = std::vector<double>(x_len, 0.0);
		myCatheter = Catheter(outside, bladder, inside, skin_concentration, bag_concentration);
		break;
	default:
		myCatheter = Catheter(skin_concentration, bag_concentration, x_len); 
	}

	BasicParameters myParam = BasicParameters(diffusivity, growth_rate1, // Parameters needed to solve the PDE problem
		carrying_capacity1, growth_rate2, carrying_capacity2, growth_rate3,
		carrying_capacity3, urine_rate, catheter_radius, stickiness, sump_volume, catheter_length);

	PDE myPDE = PDE(&myParam, &myCatheter, simulation_length, dt, print_step);
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
		results_file <<"simulation length,time step,print step,diffusivity,outside growth rate,"
			<< "outside carrying capacity,bladder growth rate,bladder carrying capacity," << 
			"inside growth rate,inside carrying capacity,catheter radius,stickiness,sump volume,"
			<< "catheter length,initial condition,num of x steps,skin concentration," <<
			"drainage bag concentration" << "\n";
		results_file << simulation_length << "," << dt << "," << print_step << "," << diffusivity <<
			"," << growth_rate1 << "," << carrying_capacity1 << "," << growth_rate2 << "," <<
			carrying_capacity2 << "," << growth_rate3 << "," << carrying_capacity3 << ","  <<
			catheter_radius << "," << stickiness << "," << sump_volume << "," << catheter_length <<
			"," << initial_condition << "," << x_len << "," << skin_concentration << "," << 
			bag_concentration << "\n\n";
		myPDE.solve(results_file);
	}

	results_file.close();

}
