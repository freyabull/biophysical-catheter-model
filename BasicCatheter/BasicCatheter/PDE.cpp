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
	std::cout << "Skin concentration is: " << data->skin_concentration << "and bag concentration is: " << data->bag_concentration << std::endl;
	// Check if there is external contamination from the skin
	if (data->skin_concentration < 0)
	{
		// Check if there is external contamination from drainage bag
		if (data->bag_concentration < 0)
		{
			std::cout << "No contaminants" << std::endl;
			// No contaminants
			for (int i = 0; i < N; ++i)
			{
				if (i%print_step == 0)
				{
					record(i);
				}
				// Deal with boundary condition
				step_c();
				step_dc();
				// Solve for current timestep
				step();
			}
	    }
		else
		{
			std::cout << "Only drainage bag contaminates" << std::endl;
			// Drainage bag contaminates
			for (int i = 0; i < N; ++i)
			{
				if (i%print_step == 0)
				{
					record(i);
				}
				// Deal with boundary condition
				step_c();
				step_de();
				// Solve for current timestep
				step();
			}
		}
	}
	else 
	{
		// Check if there is external contamination from the drainage bag
		if (data->bag_concentration < 0)
		{
			std::cout << "Only skin contaminates" << std::endl;
			// Skin contamination
			for (int i = 0; i < N; ++i)
			{
				if (i%print_step == 0)
				{
					record(i);
				}
				// Deal with boundary condition
				step_e();
				step_dc();
				// Solve for current timestep
				step();
			}
		}
		else
		{
			std::cout << "Both skin and drainage contaminates" << std::endl;
			// Both bag and skin contamination
			for (int i = 0; i < N; ++i)
			{
				if (i%print_step == 0)
				{
					record(i);
				}
				// Deal with boundary condition
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
				step_c();
				step_de();
				// Solve for current timestep
				step();
			}
		}
	}
	else
	{
		// Check if there is external contamination from the drainage bag
		if (data->bag_concentration < 0)
		{
			// Skin contamination
			for (int i = 0; i < N; ++i)
			{
				if (i%print_step == 0)
				{
					record(i, file);
				}
				// Deal with boundary condition
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
	data->outside[x_len-1] = co2 * data->old_outside[x_len-2] + co3*data->bladder + coi*data->old_inside[0]
		+ data->old_outside[x_len-1] * (co4 - co5 * data->old_outside[x_len-1]);
	// Fisher wave equation for outside of catheter
	for (int j = 1; j < x_len-1; ++j) 
	{
		data->outside[j] = co1 * (data->old_outside[j+1]
			+ data->old_outside[j-1]) + data->old_outside[j] 
			* (co4 - co5 * data->old_outside[j]);
	}

	
	// Boundary conditions. Diffusion across from bladder to inside
	data->inside[0] = ci2 * data->old_inside[1] + ci3 * data->bladder + coi*data->old_outside[x_len-1]
		+ data->old_inside[0] * (ci4 - ci5 * data->old_inside[0]) + ci7;
	// Fisher wave equation with source term for inside of catheter
	for (int k = 1; k < x_len-1; ++k) 
	{
		data->inside[k] = ci1 * (data->old_inside[k + 1]
			+ data->old_inside[k - 1]) + data->old_inside[k] * (ci4 - ci5
				* data->old_inside[k]) + ci7;
	}

	// Bladder is a well-mixed volume with diffusion to catheter
	data->bladder = cb1 * data->old_outside[x_len - 1] + cb2 *
		data->old_inside[0] + data->bladder * (cb3 - cb4 * data->bladder);

	// Update old data
	data->update();
	ci7 = ci6 * data->bladder; // Update constant governing intraluminal deposition to include current concentration
}

void PDE::step_c()
{
	// No external contamination boundary condition: no flux.
	data->outside[0] = co1 * (2.0 * data->old_outside[1])
		+ data->old_outside[0] * (co4 - co5 * data->old_outside[0]);
}

void PDE::step_e()
{
	// External contamination from skin: fixed boundary.
	data->outside[0] = co1 * (data->old_outside[1] + data->skin_concentration)
		+ data->old_outside[0] * (co4 - co5 * data->old_outside[0]);
}

void PDE::step_dc()
{
	//  No external contamination boundary condition: no flux.
	data->inside[x_len - 1] = ci1 * (2.0 * data->old_inside[x_len - 2])
		+ data->old_inside[x_len - 1] * (ci4 - ci5 * data->old_inside[x_len - 1])
		+ ci7;
}

void PDE::step_de()
{
	//  External contamination boundary condition: fixed boundary.
	data->inside[x_len - 1] = ci1 * (data->old_inside[x_len - 2] 
		+ data->bag_concentration) + data->old_inside[x_len - 1] * (ci4 - ci5 
			* data->old_inside[x_len - 1]) + ci7;
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
	x_len = data->x_len;
	N = int(time / dt);
	print_step = int(print_interval / dt);
	dx = param->catheter_length / (x_len-1);
	co1 = param->diffusivity * dt / (dx*dx);
	co2 = param->diffusivity * dt * (2 - 0.5*param->coupling_o - 0.5*param->coupling_oi) / (dx*dx);
	co3 = param->diffusivity * dt * param->coupling_o / (dx*dx);
	coi = param->diffusivity * dt * param->coupling_oi / (dx*dx);
	co4 = 1.0 - 2.0 * co1 + param->growth_rate1 * dt;
	co5 = param->growth_rate1 * dt / param->carrying_capacity1;
	cb1 = param->diffusivity * dt * param->coupling_o/ (param->sump_volume);
	cb2 = param->diffusivity * dt * param->coupling_i/ (param->sump_volume);
	cb3 = 1.0 + param->growth_rate2 * dt - (param->coupling_o+param->coupling_i) * param->diffusivity * dt / param->sump_volume;
	cb4 = param->growth_rate2 * dt / param->carrying_capacity2;
	ci1 = param->diffusivity * dt / (dx*dx);
	ci2 = param->diffusivity * dt * (2 - 0.5*param->coupling_i - 0.5*param->coupling_oi) / (dx*dx);
	ci3 = param->diffusivity * dt * param->coupling_i / (dx*dx);
	ci4 = 1.0 - 2.0 * ci1 + param->growth_rate3*dt;
	ci5 = param->growth_rate3 * dt / param->carrying_capacity3;
	ci6 = param->stickiness * 2 * 3.14 * shell_thickness * param->catheter_radius * dt;
	ci7 = ci6 * data->bladder;
}