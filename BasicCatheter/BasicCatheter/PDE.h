#pragma once
#include "BasicParameters.h"
#include "Catheter.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include "Integrate.h"
#include <vector>

class PDE // Methods for solving the pde problem
{
public:
	BasicParameters* param; // Parameters needed to solve the PDE problem
	Catheter* data; // Current state of catheter
	int x_len; // Number of x points 
	int r_len; // Number of r points
	double time; // Time of total simulation
	double dt; // Time-step
	double print_interval; // Time interval at which to record data
	double dx; // x step
	double dr; // r step

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
	double R_square; // Internal catheter radius squared
	double dr_square; // r step squared
	double o1; // Constant for outside of catheter: D dt/ dx^2
	double o2; // Constant for outside of catheter: 1 - 2 D dt / dx^2 + r dt
	double o3; // Constant for outside of catheter: r dt / kappa
	double o4; // Constant for outside of catheter: k_a V_c dt / (2 pi R_e dx)
	double o5; // Constant for outside of catheter: 1 - 2 D dt / dx^2 + r dt - k_d S_c dt / (2 pi R_e dx)
	double coi; // Constant for diffusion from outside to inside: D dt / dx^2
	double b1; // Constant for bladder: 1 + r dt - k_D dt / V - k_a V_c dt / V
	double b2; // Constant for bladder: r dt / kappa
	double b3; // Constant for bladder: k_d S_c dt / V
	double i1; // Constant for inside of catheter: D dt / dx^2
	double i2; // Constant for inside of catheter: 1 - 2 D dt / dx^2 + r dt
	double i3; // Constant for inside of catheter: r dt / kappa
	double f1; // Constant for inside flow : D dt / dr^2
	double f2; // Constant for inside flow : 2 lambda dt / (pi R^4 dx)
	double f3; // Constant for inside flow : 1 - 2 D dt / dr^2
	double f7b; // Constant for inside flow : 2 lambda dt / (pi R^2 dx)
	double f8b; // Constant for inside flow : 1 - 2 D dt / dr^2 - 2 lambda dt / (pi R^2 dx)
	double f9; // Constant for inside flow : D dt / dr
	double of1; // Constant for calculation of outflow density : 4 dr / R^2
	double of2; // Constant for calculation of outflow density : 4 dr^3 / R^4
	// Solve for current time-step using explicit methods
	void step(); 
	// Case for no external contamination (skin_concentration < 0)
	void step_bc_skin_clean(); 
	// Case for external contamination (skin_concentration > 0)
	void step_bc_skin_contamination(); 
	// Case for no drainage contamination (bag_concentration < 0)
	void step_bc_drainage_clean();
	// Case for drainage contamination (bag_concentration > 0)
	void step_bc_drainage_contamination();
	// Find the average density in the flow out the bottom of the catheter
	void find_outflow_density();
	// Output results for current time-step to console
	void record(int current_step); 
	// Output results for current time-step to given file
	void record(int current_step, std::ofstream &file); 
	// Initialize variables
	void initialize(); 
	// Check for numerical stability
	void stability_check();
};

