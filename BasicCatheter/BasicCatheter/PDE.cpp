// Methods for solving the PDE problem

#include "pch.h"
#include "PDE.h"
#include <iostream>


PDE::PDE()
{
}

PDE::PDE(BasicParameters* t_param, Catheter* t_data, double t_time,
	     double t_dt, double t_print_interval) :
	param(t_param), data(t_data), time(t_time), dt(t_dt), 
	print_interval(t_print_interval)
{
	initialize();
}

PDE::~PDE()
{
}

void PDE::solve()
{
	// Check if there is external contamination from the skin
	if (data->skin_concentration < 0) {
		for (int i = 0; i < N; ++i) {
			if (i%print_step == 0) {
				std::cout << "Time is " << i * dt / 3600 << "hrs \n";
				record();
			}
			// Deal with boundary condition
			step_c();
			// Solve for current timestep
			step();
		}
	}
	else {
		for (int i = 0; i < N; ++i) {
			if (i%print_step == 0) {
				std::cout << "Time is " << i * dt / 3600 << "hrs \n";
				record();
			}
			// Deal with boundary condition
			step_e();
			// Solve for current timestep
			step();
		}
	}
}

void PDE::step()
{
	// Outside
	// top of catheter
	data->outside[x1_len-1] = c1 * (data->bladder + data->old_outside[x1_len-2])
		+ data->old_outside[x1_len-1] * (c2 - c3 * data->old_outside[x1_len-1]);
	// main steps
	for (int j = 1; j < x1_len-1; ++j) {
		data->outside[j] = c1 * (data->old_outside[j+1] + data->old_outside[j-1])
			+ data->old_outside[j] * (c2 - c3 * data->old_outside[j]);
	}
	// Bladder
	data->bladder = c4 * data->old_outside[x1_len-1] + data->bladder * 
		            (c5 - c6 * data->bladder);
	// Inside
	for (int k = 0; k < x2_len; ++k) {

	}
	// update old data
	data->update();
}

void PDE::step_c()
{
	// no external contamination
	data->outside[0] = c1 * (2 * data->old_outside[1])
		+ data->old_outside[0] * (c2 - c3 * data->old_outside[0]);
}

void PDE::step_e()
{
	// external contamination from skin
	data->outside[0] = c1 * (data->old_outside[1] + data->skin_concentration)
		+ data->old_outside[0] * (c2 - c3 * data->old_outside[0]);
}

void PDE::record()
{ /*
	std::cout << "Bottom concentration is " << data->outside[0] << "\n";
	std::cout << "Top concentration is " << data->outside[x1_len-1] << "\n";
	std::cout << "Bladder concentration is " << data->bladder << "\n";
*/}

void PDE::initialize()
{
	x1_len = data->outside.size();
	x2_len = data->inside.size();
	N = int(time / dt);
	print_step = int(print_interval / dt);
	dx1 = param->catheter_length / x1_len;
	dx2 = param->catheter_length / x2_len;
	c1 = param->diffusivity * dt / (dx1*dx1);
	c2 = 1 - 2 * c1 + param->growth_rate1 * dt;
	c3 = param->growth_rate1 * dt / param->carrying_capacity1;
	c4 = 2 * param->diffusivity * dt / (dx1*dx1*param->sump_volume);
	c5 = param->growth_rate2 * dt + 1;
	c6 = param->growth_rate2 * dt / param->carrying_capacity2;
}