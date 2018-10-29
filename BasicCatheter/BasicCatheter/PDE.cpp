// Methods for solving the PDE problem

#include "pch.h"
#include "PDE.h"



PDE::PDE()
{
}

PDE::PDE(BasicParameters* t_param, Catheter* t_data,
	double t_time, double t_dt, double t_print_interval) :
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
	if (data->skin_concentration < 0) 
	{
		for (int i = 0; i < N; ++i) 
		{
			if (i%print_step == 0) 
			{
				record(i);
			}
			// Deal with boundary condition
			step_c();
			// Solve for current timestep
			step();
		}
	}
	else 
	{
		for (int i = 0; i < N; ++i) 
		{
			if (i%print_step == 0) 
			{
				record(i);
			}
			// Deal with boundary condition
			step_e();
			// Solve for current timestep
			step();
		}
	}
}


void PDE::solve(std::ofstream &file)
{
	// Check if there is external contamination from the skin
	if (data->skin_concentration < 0)
	{
		for (int i = 0; i < N; ++i)
		{
			if (i%print_step == 0)
			{
				record(i,file);
			}
			// Deal with boundary condition
			step_c();
			// Solve for current timestep
			step();
		}
	}
	else
	{
		for (int i = 0; i < N; ++i)
		{
			if (i%print_step == 0)
			{
				record(i,file);
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
	// Boundary condition for top of catheter (diffusion into bladder)
	data->outside[x1_len-1] = c1 * (data->bladder 
		+ data->old_outside[x1_len-2]) + data->old_outside[x1_len-1]
		* (c2 - c3 * data->old_outside[x1_len-1]);
	// Fisher wave equation for outside of catheter
	for (int j = 1; j < x1_len-1; ++j) 
	{
		data->outside[j] = c1 * (data->old_outside[j+1]
			+ data->old_outside[j-1]) + data->old_outside[j] 
			* (c2 - c3 * data->old_outside[j]);
	}

	// Update c10 to include current concentration
	double c11 = c10 * data->bladder;
	// Boundary conditions. No flux conditions at both ends
	data->inside[0] = c7 * (2.0 * data->old_inside[1]) + data->old_inside[0]
		* (c8 - c9 * data->old_inside[0]) + c11 * pow(1.0, 0.5);
	data->inside[x2_len-1] = c7 * (2.0 * data->old_inside[x2_len-2]) 
		+ data->old_inside[x2_len-1] * (c8 - c9 * data->old_inside[x2_len-1])
		+ c11 * pow(1.0 / x2_len, 0.5);
	// Fisher wave equation with source term for inside of catheter
	for (int k = 1; k < x2_len-1; ++k) 
	{
		data->inside[k] = c7 * (data->old_inside[k + 1]
			+ data->old_inside[k - 1]) + data->old_inside[k] * (c8 - c9
				* data->old_inside[k]) + c11 * 1.0/pow(double(k)+1.0, 0.5);
	}

	// Bladder is a well-mixed volume with diffusion to catheter
	data->bladder = c4 * data->old_outside[x1_len - 1] + c4a * 
		data->old_inside[0] + data->bladder * (c5 - c6 * data->bladder);

	// Update old data
	data->update();
}

void PDE::step_c()
{
	// No external contamination boundary condition: no flux.
	data->outside[0] = c1 * (2.0 * data->old_outside[1])
		+ data->old_outside[0] * (c2 - c3 * data->old_outside[0]);
}

void PDE::step_e()
{
	// External contamination from skin: fixed boundary.
	data->outside[0] = c1 * (data->old_outside[1] + data->skin_concentration)
		+ data->old_outside[0] * (c2 - c3 * data->old_outside[0]);
}

void PDE::record(int current_step)
{ 
		std::cout << "Time is " << current_step * dt / 3600 << "hrs \n";
		std::cout << "Bottom conc is " << data->outside[0] << "\n";
		std::cout << "Top conc is " << data->outside[x1_len-1] << "\n";
		std::cout << "Bladder conc is " << data->bladder << "\n";
		std::cout << "Top inside conc is " << data->inside[0] << "\n";
		std::cout << "Bottom inside conc is " << data->inside[x2_len-1] << "\n";
}

void PDE::record(int current_step, std::ofstream &file)
{
	file << "\n Time is " << current_step * dt / 3600 << "hrs \n";
	file << "Outside of catheter starting at base. \n";
	for (int i = 0; i < x1_len; ++i) {
		file << data->outside[i] << ",";
	}
	file << "\n Bladder \n";
	file << data->bladder << "\n";
	file << "Inside of catheter starting at top. \n";
	for (int i = 0; i < x1_len; ++i) {
		file << data->inside[i] << ",";
	}
	file << "\n";
}

void PDE::initialize()
{
	x1_len = data->outside.size();
	x2_len = data->inside.size();
	N = int(time / dt);
	print_step = int(print_interval / dt);
	dx1 = param->catheter_length / x1_len;
	dx2 = param->catheter_length / x2_len;
	double num_drop = 60.0 * param->droplet_size / param->urine_rate; // Number of droplets per second
	c1 = param->diffusivity * dt / (dx1*dx1);
	c2 = 1.0 - 2.0 * c1 + param->growth_rate1 * dt;
	c3 = param->growth_rate1 * dt / param->carrying_capacity1;
	c4 = param->diffusivity * dt / (dx1*dx1*param->sump_volume);
	c4a = param->diffusivity * dt / (dx2*dx2*param->sump_volume);
	c5 = 1.0 - c4 - c4a - param->growth_rate2 * dt;
	c6 = param->growth_rate2 * dt / param->carrying_capacity2;
	c7 = param->diffusivity *dt / (dx2*dx2);
	c8 = 1.0 - 2.0 * c7 + param->growth_rate3*dt;
	c9 = param->growth_rate3 *dt / param->carrying_capacity3;
	c10 = pow(3.14, 5.0 / 3.0) * pow(0.75*param->droplet_size, 2.0 / 3.0)
		 * pow(dx2 / (2.0 * 9.81e3), 0.5) * param->stickiness * dt/num_drop;
}