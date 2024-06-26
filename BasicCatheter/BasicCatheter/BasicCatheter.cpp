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
#include <algorithm>
#include "Catheter.h"
#include "BasicParameters.h"
#include "PDE.h"
#include <time.h>

int main(int argc, char **argv)
{
	clock_t start, end;
	start = clock();
    // Define variables
	double diffusivity = 1e-4; // diffusivity of bacteria in urine, in mm^2/s
	double surface_diffusivity = 1e-7; // diffusivity of bacteria on catheter, in mm^2/s
	double growth_rate1 = std::log(2.0) / 3600.0; // Growth rate of bacteria on outside of catheter, per s
	double carrying_capacity1 = 1e7; // Carrying capacity of bacteria on outside of catheter, per mm^2
	double growth_rate2 = std::log(2.0) / 1800.0; // Growth rate of bacteria in bladder, per s
	double carrying_capacity2 = 1e6; // Carrying capacity of bacteria in bladder, per mm^3 
	double growth_rate3 = std::log(2.0) / 3600.0; // Growth rate of bacteria on inside of catheter, per s
	double carrying_capacity3 = 1e7; // Carrying capacity of bacteria on inside of catheter, per mm^2
	double urine_rate = 5.0;//1.0e3 / 60.0; // Rate of generation of urine by kidneys, in mm^3/s 
	double catheter_radius = 1.0; // Radius of catheter in mm
	double stickiness = 1.0; // Probability of a bacterium sticking to surface if it comes into contact
	double sump_volume = 5.0e4; // Volume of residual urine in bladder, in mm^3
	double catheter_length = 40.0; // Length of catheter, in mm
	double skin_concentration = 1e2; // Concentration of bacteria on skin (ie boundary condition)
	double bag_concentration = -0.1; // Concentration of bacteria in drainage bag (ie boundary condition)
	int initial_condition = 0; // Type of initial conditions
	int x_len = catheter_length * 170; // Number of x points
	int r_len = 11; // Number of r points
	const char* file_name = "results.csv"; // File for output to be written to
	double initial_skin_conc = 1e2; // Concentration of bacteria on the skin at time of catheter insertion
	double lambda = 1.0/catheter_length; // Scale factor for exponential distributions
	int simulation_length = 86400 * 200; // Timeframe of simulation, in s
	double dt = 60; // Time step
	int print_interval = 3600; // Time interval at which to output data (s)
	int print_num_steps = 101; // Number of x steps at which to output data
	double catheter_external_radius = catheter_radius;// +1.0; // External catheter radius in mm
	double attachment_rate = 4*3.14*diffusivity*1e-3; // Rate at which a bacterium in contact sticks (s^-1) Smoluchowski - 4 pi D sigma
	double detachment_rate = growth_rate1; // Rate at which a bacterium in contact detaches (s^-1) say all new cells formed at boundary detach
	double viscosity = 0.83; // Viscosity of urine (mm^2 s^-1)


	// Set initial conditions
	Catheter myCatheter; // Current state of catheter
	std::vector<double> outside; // Concentration profile of outside of catheter
	double bladder; // Concentration within bladder
	std::vector<double> inside; // Concentration profile of inside of catheter

	switch (initial_condition) 
	{
	case 1:
		// Clean catheter
		myCatheter = Catheter(skin_concentration, bag_concentration, x_len, r_len); 
		break;
	case 2:
		// Uniform distribution
		outside = std::vector<double>(x_len, initial_skin_conc/(x_len-1));
		bladder = 0.0;
		inside = std::vector<double>(x_len, 0.0);
		myCatheter = Catheter(outside, bladder, inside, skin_concentration, bag_concentration, r_len);
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
		myCatheter = Catheter(outside, bladder, inside, skin_concentration, bag_concentration, r_len);
		break;
	case 4:
		// Infected bladder
		outside = std::vector<double>(x_len, 0.0);
		bladder = 1e2;
		inside = std::vector<double>(x_len, 0.0);
		myCatheter = Catheter(outside, bladder, inside, skin_concentration, bag_concentration, r_len);
		break;
	default:
		myCatheter = Catheter(skin_concentration, bag_concentration, x_len, r_len); 
	}
	BasicParameters myParam = BasicParameters(diffusivity, surface_diffusivity, growth_rate1, // Parameters needed to solve the PDE problem
		carrying_capacity1, growth_rate2, carrying_capacity2, growth_rate3,
		carrying_capacity3, urine_rate, catheter_radius, catheter_external_radius, stickiness, 
		sump_volume, catheter_length, viscosity, attachment_rate, detachment_rate );

	PDE myPDE = PDE(&myParam, &myCatheter, simulation_length, dt, print_interval, print_num_steps);
	std::ofstream results_file; // File for output to be written to 
	results_file.open(file_name, std::ios::trunc);
	if (!results_file.is_open())
	{
		std::cout << "Unable to open file " << file_name << std::endl;
		std::cout << "Writing direct to console instead" << std::endl;
		myPDE.solve_light();
	}
	else
	{ 
		std::cout << "Successfully opened " << file_name << std::endl;
		results_file <<"simulation length,time step,print interval,print num steps,diffusivity,surface diffusivity,outside growth rate,"
			<< "outside carrying capacity,bladder growth rate,bladder carrying capacity," << 
			"inside growth rate,inside carrying capacity,urine rate,catheter radius,external catheter radius,stickiness,sump volume,"
			<< "catheter length,initial condition,num of x steps,skin concentration," <<
			"drainage bag concentration,viscosity,attachment rate,detachment rate" << "\n";
		results_file << simulation_length << "," << dt << "," << print_interval << "," <<print_num_steps<<","<< diffusivity << "," << surface_diffusivity <<
			"," << growth_rate1 << "," << carrying_capacity1 << "," << growth_rate2 << "," <<
			carrying_capacity2 << "," << growth_rate3 << "," << carrying_capacity3 << ","  << urine_rate << "," <<
			catheter_radius << "," << catheter_external_radius << "," << stickiness << "," << sump_volume << "," << catheter_length <<
			"," << initial_condition << "," << x_len << "," << skin_concentration << "," << 
			bag_concentration << "," << viscosity << "," << attachment_rate << "," << detachment_rate << "\n\n";
		clock_t pre_solve_time = clock();
		std::cout << "Time elapsed in set up " << double(pre_solve_time - start) / CLOCKS_PER_SEC << std::endl;
		myPDE.solve_light(results_file);
	}

	results_file.close();
	end = clock();

	double duration_sec = double(end - start) / CLOCKS_PER_SEC;
	std::cout << "Time elapsed " << duration_sec << std::endl;

}