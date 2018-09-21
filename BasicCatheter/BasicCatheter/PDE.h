#pragma once
#include "BasicParameters.h"
#include "Catheter.h"
class PDE // Methods for solving the pde problem
{
public:
	BasicParameters* param;
	Catheter* data;
	int x1_len; // Number of x steps on outside of catheter
	int x2_len; // Number of x steps on inside of catheter
	double time; // Time of total simulation
	double dt; // Time-step
	double print_interval; // Time interval at which to record data
	double dx1; // x step for outside of catheter
	double dx2; // x step for inside of catheter

	PDE();
	PDE(BasicParameters* param, Catheter* catheter, double time = 3600,
		double dt = 0.1, double print_interval = 600);
	~PDE();

	void solve(); // Solve the model for the parameters given

private:
	int N; // Number of time steps
	int print_step; // Time-step at which to output results
	double c1; // Constant for outside of catheter
	double c2; // Constant for outside of catheter
	double c3; // Constant for outside of catheter
	double c4; // Constant for bladder
	double c5; // Constant for bladder
	double c6; // Constant for bladder
	void step(); // Solve for current time-step
	void step_c(); // Case for no external contamination (skin_concentration < 0)
	void step_e(); // Case for external contamination (skin_concentration > 0)
	void record(); // Output results for current time-step
	void initialize(); // Initialize variables
};

