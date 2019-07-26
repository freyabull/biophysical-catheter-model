#pragma once
#include "BasicParameters.h"
#include "Catheter.h"
#include <iostream>
#include <cmath>
#include <fstream>

class PDE // Methods for solving the pde problem
{
public:
	BasicParameters* param; // Parameters needed to solve the PDE problem
	Catheter* data; // Current state of catheter
	int x_len; // Number of x points 
	double time; // Time of total simulation
	double dt; // Time-step
	double print_interval; // Time interval at which to record data
	double dx; // x step

	PDE();
	PDE(BasicParameters* param, Catheter* catheter,
		double time = 3600, double dt = 0.1, double print_interval = 600);
	~PDE();

	// Solve the model for the parameters set. Writes output to console
	void solve();
	// Solve the model for the parameters set. Writes output to file
	void solve(std::ofstream &file); 

private:
	int N; // Number of time steps
	int print_step; // Time-step at which to output results
	double co1; // Constant for outside of catheter
	double co2; // Constant for outside of catheter
	double co3; // Constant for outside of catheter
	double coi; // Constant for diffusion between inside and outside
	double co4; // Constant for outside of catheter
	double co5; // Constant for outside of catheter
	double cb1; // Constant for bladder
	double cb2; // Constant for bladder
	double cb3; // Constant for bladder
	double cb4; // Constant for bladder
	double ci1; // Constant for inside of catheter
	double ci2; // Constant for inside of catheter
	double ci3; // Constant for inside of catheter
	double ci4; // Constant for inside of catheter
	double ci5; // Constant for inside of catheter
	double ci6; // Coefficient of intraluminal deposition rate
	double ci7; // Intraluminal deposition rate
	// Solve for current time-step
	void step(); 
	// Case for no external contamination (skin_concentration < 0)
	void step_c(); 
	// Case for external contamination (skin_concentration > 0)
	void step_e(); 
	// Case for no drainage contamination (bag_concentration < 0)
	void step_dc();
	// Case for drainage contamination (bag_concentration > 0)
	void step_de();
	// Output results for current time-step to console
	void record(int current_step); 
	// Output results for current time-step to given file
	void record(int current_step, std::ofstream &file); 
	// Initialize variables
	void initialize(); 
};

