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
		// Check if there is external contamination from drainage bag
		if (data->bag_concentration < 0)
		{
			// No contaminants
			for (int i = 0; i < N; ++i)
			{
				if (i%print_step == 0)
				{
					record(i);
				}
				// Deal with boundary condition
				c11 = c10 * data->bladder; // Update c10 to include current concentration
				step_c();
				step_dc();
				// Solve for current timestep
				step();
			}
	    }
		else
		{
			// Drainage bag contaminates
			for (int i = 0; i < N; ++i)
			{
				if (i%print_step == 0)
				{
					record(i);
				}
				// Deal with boundary condition
				c11 = c10 * data->bladder; // Update c10 to include current concentration
				step_c();
				step_de();
				// Solve for current timestep
				step();
			}
		}
	}
	else 
	{
		// Check if there is external contamination from the skin
		if (data->skin_concentration < 0)
		{
			// Skin contamination
			for (int i = 0; i < N; ++i)
			{
				if (i%print_step == 0)
				{
					record(i);
				}
				// Deal with boundary condition
				c11 = c10 * data->bladder; // Update c10 to include current concentration
				step_e();
				step_dc();
				// Solve for current timestep
				step();
			}
		}
		else
		{
			// Both bag and skin contamination
			for (int i = 0; i < N; ++i)
			{
				if (i%print_step == 0)
				{
					record(i);
				}
				// Deal with boundary condition
				c11 = c10 * data->bladder; // Update c10 to include current concentration
				step_e();
				step_de();
				// Solve for current timestep
				step();
			}
		}
	}
}


void PDE::solve(std::ofstream &file)
{
	// Check if there is external contamination from the skin
	if (data->skin_concentration < 0)
	{
		// Check if there is external contamination from drainage bag
		if (data->bag_concentration < 0)
		{
			// No contaminants
			for (int i = 0; i < N; ++i)
			{
				if (i%print_step == 0)
				{
					record(i, file);
				}
				// Deal with boundary condition
				c11 = c10 * data->bladder; // Update c10 to include current concentration
				step_c();
				step_dc();
				// Solve for current timestep
				step();
			}
		}
		else
		{
			// Drainage bag contaminates
			for (int i = 0; i < N; ++i)
			{
				if (i%print_step == 0)
				{
					record(i, file);
				}
				// Deal with boundary condition
				c11 = c10 * data->bladder; // Update c10 to include current concentration
				step_c();
				step_de();
				// Solve for current timestep
				step();
			}
		}
	}
	else
	{
		// Check if there is external contamination from the skin
		if (data->skin_concentration < 0)
		{
			// Skin contamination
			for (int i = 0; i < N; ++i)
			{
				if (i%print_step == 0)
				{
					record(i, file);
				}
				// Deal with boundary condition
				c11 = c10 * data->bladder; // Update c10 to include current concentration
				step_e();
				step_dc();
				// Solve for current timestep
				step();
			}
		}
		else
		{
			// Both bag and skin contamination
			for (int i = 0; i < N; ++i)
			{
				if (i%print_step == 0)
				{
					record(i, file);
				}
				// Deal with boundary condition
				c11 = c10 * data->bladder; // Update c10 to include current concentration
				step_e();
				step_de();
				// Solve for current timestep
				step();
			}
		}
	}
}

void PDE::step()
{
	// Boundary condition for top of catheter (diffusion into bladder)
	data->outside[x_len-1] = c1 * (data->bladder 
		+ data->old_outside[x_len-2]) + data->old_outside[x_len-1]
		* (c2 - c3 * data->old_outside[x_len-1]);
	// Fisher wave equation for outside of catheter
	for (int j = 1; j < x_len-1; ++j) 
	{
		data->outside[j] = c1 * (data->old_outside[j+1]
			+ data->old_outside[j-1]) + data->old_outside[j] 
			* (c2 - c3 * data->old_outside[j]);
	}

	
	// Boundary conditions. Diffusion across from bladder
	data->inside[0] = c7 * (data->bladder + data->old_inside[1]) + data->old_inside[0]
		* (c8 - c9 * data->old_inside[0]) + c11;
	// Fisher wave equation with source term for inside of catheter
	for (int k = 1; k < x_len-1; ++k) 
	{
		data->inside[k] = c7 * (data->old_inside[k + 1]
			+ data->old_inside[k - 1]) + data->old_inside[k] * (c8 - c9
				* data->old_inside[k]) + c11;
	}

	// Bladder is a well-mixed volume with diffusion to catheter
	data->bladder = c4 * (data->old_outside[x_len - 1] +  
		data->old_inside[0]) + data->bladder * (c5 - c6 * data->bladder);

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

void PDE::step_dc()
{
	//  No external contamination boundary condition: no flux.
	data->inside[x_len - 1] = c7 * (2.0 * data->old_inside[x_len - 2])
		+ data->old_inside[x_len - 1] * (c8 - c9 * data->old_inside[x_len - 1])
		+ c11;
}

void PDE::step_de()
{
	//  External contamination boundary condition: fixed boundary.
	data->inside[x_len - 1] = c7 * (data->old_inside[x_len - 2] 
		+ data->bag_concentration) + data->old_inside[x_len - 1] * (c8 - c9 
			* data->old_inside[x_len - 1]) + c11;
}

void PDE::record(int current_step)
{ 
	std::cout << "Time is " << double(current_step) * dt / 3600.0 << "hrs \n";
	std::cout << "Bottom conc is " << data->outside[0] << "\n";
	std::cout << "Top conc is " << data->outside[x_len-1] << "\n";
	std::cout << "Bladder conc is " << data->bladder << "\n";
	std::cout << "Top inside conc is " << data->inside[0] << "\n";
	std::cout << "Bottom inside conc is " << data->inside[x_len-1] << "\n\n";
}

void PDE::record(int current_step, std::ofstream &file)
{
	double out_time = double(current_step) * dt / 3600.0;
	file << "o" << out_time;
	for (int i = 0; i < x_len; ++i) {
		file << "," << data->outside[i];
	}
	file << "\nb" << out_time;
	file << "," << data->bladder << "\n";
	file << "i" << out_time ;
	for (int i = 0; i < x_len; ++i) {
		file << "," << data->inside[i];
	}
	file << "\n";
}

void PDE::initialize()
{
	double shell_thickness = 5e-3; // Width of contact area
	x_len = data->outside.size();
	N = int(time / dt);
	print_step = int(print_interval / dt);
	dx = param->catheter_length / x_len;
	c1 = param->diffusivity * dt / (dx*dx);
	c2 = 1.0 - 2.0 * c1 + param->growth_rate1 * dt;
	c3 = param->growth_rate1 * dt / param->carrying_capacity1;
	c4 = param->diffusivity * dt / (param->sump_volume);
	c5 = 1.0 + param->growth_rate2 * dt;
	c6 = param->growth_rate2 * dt / param->carrying_capacity2;
	c7 = param->diffusivity *dt / (dx*dx);
	c8 = 1.0 - 2.0 * c7 + param->growth_rate3*dt;
	c9 = param->growth_rate3 * dt / param->carrying_capacity3;
	c10 = param->stickiness * 2 * 3.14 * shell_thickness * param->catheter_radius * dt;
	c11 = c10 * data->bladder;
}