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
	std::string file_name; // File for output to be written to
	double catheter_length; // Length of catheter, in mm
	double sump_volume; // Volume of residual urine in bladder, in mm^3
	double urine_rate; // Rate of generation of urine by kidneys, in mm^3/s 
	double catheter_radius; // Radius of catheter in mm
	double surface_diffusivity; // diffusivity of bacteria on catheter, in mm^2/s
	double diffusivity; // diffusivity of bacteria in urine, in mm^2/s
	double growth_rate1; // Growth rate of bacteria on outside of catheter, per s
	double growth_rate2; // Growth rate of bacteria in bladder, per s
	double carrying_capacity1; // Carrying capacity of bacteria on outside of catheter, per mm^2
	double carrying_capacity2; // Carrying capacity of bacteria in bladder, per mm^3 

	// Open the input file
	const char* param_file_name = "params.txt";
	std::ifstream param_file(param_file_name);
	if (!param_file.is_open()) {
		std::cout << "Unable to open parameter file " << param_file_name << " proceeding with preset parameters" << std::endl;
		// Set variables from preset values
		file_name = "default_results.csv"; // File for output to be written to
		catheter_length = 40.0; // Length of catheter, in mm
		sump_volume = 5.0e4; // Volume of residual urine in bladder, in mm^3
		urine_rate = 1.0e3 / 60.0; // Rate of generation of urine by kidneys, in mm^3/s 
		catheter_radius = 1.0; // Radius of catheter in mm
		surface_diffusivity = 1e-6; // diffusivity of bacteria on catheter, in mm^2/s
		diffusivity = 1e-4; // diffusivity of bacteria in urine, in mm^2/s
		growth_rate1 = std::log(2.0) / 3600.0; // Growth rate of bacteria on outside of catheter, per s
		growth_rate2 = std::log(2.0) / 1800.0; // Growth rate of bacteria in bladder, per s
		carrying_capacity1 = 1e7; // Carrying capacity of bacteria on outside of catheter, per mm^2
		carrying_capacity2 = 1e6; // Carrying capacity of bacteria in bladder, per mm^3 
	} else {
		// Read parameters from input file
		param_file >> file_name;
		param_file >> catheter_length;
		param_file >> sump_volume;
		param_file >> urine_rate;
		param_file >> catheter_radius;
		param_file >> surface_diffusivity;
		param_file >> diffusivity;
		param_file >> growth_rate1;
		param_file >> growth_rate2;
		param_file >> carrying_capacity1;
		param_file >> carrying_capacity2;
	}
	param_file.close();
	

	// Fixed parameters
	double growth_rate3 = growth_rate1; // Growth rate of bacteria on inside of catheter, per s
	double carrying_capacity3 = carrying_capacity1; // Carrying capacity of bacteria on inside of catheter, per mm^2
	double stickiness = 1.0; // Probability of a bacterium sticking to surface if it comes into contact
	double skin_concentration = 1e2; // Concentration of bacteria on skin (ie boundary condition)
	double bag_concentration = -0.1; // Concentration of bacteria in drainage bag (ie boundary condition)
	int initial_condition = 0; // Type of initial conditions
	double initial_skin_conc = 1e2; // Concentration of bacteria on the skin at time of catheter insertion
	double lambda = 1.0/catheter_length; // Scale factor for exponential distributions
	int r_len = 11; // Number of r points
	int target_print_num_steps = 101; // Number of x steps at which to output data
	double catheter_external_radius = catheter_radius;// +1.0; // External catheter radius in mm
	double attachment_rate = 4*3.14*diffusivity*1e-3; // Rate at which a bacterium in contact sticks (s^-1) Smoluchowski - 4 pi D sigma
	double detachment_rate = growth_rate1; // Rate at which a bacterium in contact detaches (s^-1) say all new cells formed at boundary detach
	double viscosity = 0.83; // Viscosity of urine (mm^2 s^-1)
	
	// Calculated parameters
	int x_len = catheter_length * std::max(3,2*int(std::ceil(std::sqrt(growth_rate1/surface_diffusivity)))); // Number of x points
	double dx = catheter_length / x_len;
	std::cout << "Number of x points is  " << x_len << ", dx is " << dx << std::endl;
	int simulation_length = 86400 * 20 +2*int(std::ceil(0.5*catheter_length/std::sqrt(growth_rate1*surface_diffusivity))) + 20*int(std::ceil(1/growth_rate1+1/growth_rate2)); // Timeframe of simulation
	double dt_cond1 = dx*dx / (20*surface_diffusivity); // Condition on dt for numerical stability of catheter surface
	double dt_cond2 = sump_volume / (urine_rate + attachment_rate*3.14*(dx + 5e-3)*(catheter_external_radius + 5e-3)); // Condition on dt for numerical stability of surface-bladder coupling
	double dt_cond3 = 1 / std::abs(growth_rate2 - lambda/sump_volume);
	double dt = std::min(std::min(std::min(60.0, dt_cond3), dt_cond2),dt_cond1); // Time step
	std::cout << "dt = " << dt << " and T = " << simulation_length << std::endl;
	int print_interval = std::max(3600, int(simulation_length/(24*200))); // Time interval at which to output data (s)
	std::cout << "Stability conditions on dt: Surface " << dt_cond1 << " Bladder-surface " << dt_cond2 << " Bladder " << dt_cond3 << std::endl;
	int print_num_steps = std::min(int(x_len/(int((x_len - 1) / (target_print_num_steps - 1))))+1,x_len);
	std::cout << "Target print steps " << target_print_num_steps << ", actual print steps " << print_num_steps << std::endl;
	


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
