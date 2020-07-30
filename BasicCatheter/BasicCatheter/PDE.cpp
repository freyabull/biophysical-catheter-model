// Methods for solving the PDE problem

#include "pch.h"
#include "PDE.h"

PDE::PDE(BasicParameters* t_param, Catheter* t_data,
	double t_time, double t_dt, double t_print_interval) :
	param(t_param), data(t_data), time(t_time), dt(t_dt),
	print_interval(t_print_interval)
{
	initialize();
	stability_check();
}


PDE::~PDE()
{
}

void PDE::solve()
{
	std::cout << "Skin concentration is: " << data->skin_concentration << " and bag concentration is: " << data->bag_concentration << std::endl;
	// Check if there is external contamination from the skin
	if (data->skin_concentration < 0)
	{
		// Check if there is external contamination from drainage bag
		if (data->bag_concentration < 0)
		{
			std::cout << "No contaminants" << std::endl;
			// No contaminants
			for (int i = 0; i < N; ++i) // Time index
			{
				if (i%print_step == 0)
				{
					record(i);
				}
				step_bc_skin_clean();
				step_bc_drainage_clean();
				step();
			}
	    }
		else
		{
			std::cout << "Only drainage bag contaminates" << std::endl;
			// Drainage bag contaminates only
			for (int i = 0; i < N; ++i) // Time index
			{
				if (i%print_step == 0)
				{
					record(i);
				}
				step_bc_skin_clean();
				step_bc_drainage_contamination();
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
			// Skin contamination only
			for (int i = 0; i < N; ++i) // Time index
			{
				if (i%print_step == 0)
				{
					record(i);
				}
				step_bc_skin_contamination();
				step_bc_skin_clean();
				step();
			}
		}
		else
		{
			std::cout << "Both skin and drainage contaminates" << std::endl;
			// Both bag and skin contamination
			for (int i = 0; i < N; ++i) // Time index
			{
				if (i%print_step == 0)
				{
					record(i);
				}
				step_bc_skin_contamination();
				step_bc_drainage_contamination();
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
			for (int i = 0; i < N; ++i) // Time index
			{
				if (i%print_step == 0)
				{
					record(i, file);
				}
				step_bc_skin_clean();
				step_bc_drainage_clean();
				step();
			}
		}
		else
		{
			// Drainage bag contaminates only 
			for (int i = 0; i < N; ++i) // Time index
			{
				if (i%print_step == 0)
				{
					record(i, file);
				}
				step_bc_skin_clean();
				step_bc_drainage_contamination();
				step();
			}
		}
	}
	else
	{
		// Check if there is external contamination from the drainage bag
		if (data->bag_concentration < 0)
		{
			// Skin contamination only 
			for (int i = 0; i < N; ++i) // Time index
			{
				if (i%print_step == 0)
				{
					record(i, file);
				}
				step_bc_skin_contamination();
				step_bc_drainage_clean();
				step();
			}
		}
		else
		{
			// Both bag and skin contamination
			for (int i = 0; i < N; ++i) // Time index
			{
				if (i%print_step == 0)
				{
					record(i, file);
				}
				step_bc_skin_contamination();
				step_bc_drainage_contamination();
				step();
			}
		}
	}
}

void PDE::step()
{
	for (int j = 1; j < x_len-1; ++j) // Longitudinal index
	{
		// Fisher wave equation for outside of catheter
		/* new[j] = n[j] + dt * (r * n[j] * (1-n[j]/kappa) 
					+ D / dx^2 * (n[j+1] - 2n[j] + n[j-1]) )
		*/
		data->outside[j] = o1 * (data->old_outside[j+1]
			+ data->old_outside[j-1]) + data->old_outside[j] 
			* (o2 - o3 * data->old_outside[j]);
		// Fisher wave equation for inside of catheter with source term from flow
		/* new[j] = n[j] + dt * (r * n[j] * (1-n[j]/kappa) 
					+  D / dx^2 * (n[j+1] - 2*n[j] + n[j-1]) + source)
		   where source = -D * (flow[j,r_len-1]-flow[j,r_len-2])/dr 
						= D/dr flow[j, r_len - 2]
		*/
		data->inside[j] = i1 * (data->old_inside[j + 1]
			+ data->old_inside[j - 1]) + data->old_inside[j] * (i2 - i3
				* data->old_inside[j]) + f9 * data->get_flow(j, r_len - 2);
		// Centre of inside flow: convection-diffusion equation but reflecting boundary
		data->set_flow(j, 0, 2 * f1 * data->get_flow(j, 1) 
			+ f7b * data->get_flow(j - 1, 0) + f8b * data->get_flow(j, 0));
	}
	
	for (int r = 1; r < r_len - 1; ++r) // Radial index
	{
		double f4 = 0.5 / r; // Flow constant: 1 / 2r
		double f5 = f1 * (1 + f4); // Flow constant: D dt / dr^2 * (1 + 1/2r)
		double f6 = f1 * (1 - f4); // Flow constant: D dt / dr^2 * (1 - 1/2r)
		double f7 = f2 * (R_square - std::pow(double(r),2) * dr_square); // Flow constant: 2 lambda dt / (pi R^4 dx) * (R^2 - r^2dr^2)
		double f8 = f3 - f7; // 1 - 2 D dt / dr^2 - 2 lambda dt / (pi R^4 dx) * (R^2 - r^2dr^2)
		for (int s = 1; s < x_len; ++s) // Longitudinal index
		{
			// Convection- diffusion equation for inside flow
			/* new[s,r] = f[s,r] + dt * (D/dr^2 * (f[s,r+1] - 2*f[s,r] + f[s,r-1]) 
			   + D/(2 * r * dr^2) * (f[s,r+1] - f[s,r-1]) 
			   - 2*lambda*(R^2-j^2*dr^2)/(pi*R^4*dx) * (f[s,r]-f[s-1,r]))
			*/
			double val = f5 * data->get_flow(s, r + 1)
				+ f6 * data->get_flow(s, r - 1) + f7 * data->get_flow(s - 1, r)
				+ f8 * data->get_flow(s, r);
			data->set_flow(s, r, val);
		}
		// Boundary at top of catheter (coupling to bladder)
		data->set_flow(0, r, f5 * data->get_flow(0, r + 1) 
			+ f6 * data->get_flow(0, r - 1) + f7 * data->bladder 
			+ f8 * data->get_flow(0, r));
	}
	// Boundary condition for top of catheter (diffusion into bladder)
	// Top of outside. Diffusion into inside & coupling with bladder
	data->outside[x_len - 1] = o1 * data->old_outside[x_len - 2] 
		+ o4 * data->bladder + coi * data->old_inside[0] 
		+ data->old_outside[x_len - 1] * (o5 - o3 * data->old_outside[x_len - 1]);
	
	// Boundary conditions. Diffusion across from bladder to inside
	// Top of inside. Diffusion across from outside
	data->inside[0] = i1 * data->old_inside[1] 
		+ coi * data->old_outside[x_len - 1] + data->old_inside[0] 
		* (i2 - i3 * data->old_inside[0]) + f9 * data->get_flow(0, r_len - 2);
	// Boundary for top centre of flow (coupling to bladder)
	data->set_flow(0, 0, 2 * f1 * data->get_flow(0, 1) 
		+ f7b * data->bladder + f8b * data->get_flow(0, 0));
	// Boundary for external centre of flow
	data->set_flow(x_len - 1, 0, 2 * f1 * data->get_flow(x_len - 1, 1) 
		+ f7b * data->get_flow(x_len - 2, 0) + f8b * data->get_flow(x_len - 1, 0));
	
	// Bladder is a well-mixed volume with coupling to catheter outside
	data->bladder = data->bladder * (b1 - b2 * data->bladder) 
		+ b3 * data->old_outside[x_len - 1];

	// Update old data
	data->update();
}

void PDE::step_bc_skin_clean()
{
	// No external contamination boundary condition: no flux.
	data->outside[0] = o1 * (2.0 * data->old_outside[1])
		+ data->old_outside[0] * (o2 - o3 * data->old_outside[0]);
}

void PDE::step_bc_skin_contamination()
{
	// External contamination from skin: fixed boundary.
	data->outside[0] = o1 * (data->old_outside[1] + data->skin_concentration)
		+ data->old_outside[0] * (o2 - o3 * data->old_outside[0]);
}

void PDE::step_bc_drainage_clean()
{
	//  No external contamination, boundary condition: no flux.
	data->inside[x_len - 1] = i1 * (2.0 * data->old_inside[x_len - 2])
		+ data->old_inside[x_len - 1] * (i2 - i3 * data->old_inside[x_len - 1])
		+ f9 * data->get_flow(x_len - 1, r_len - 2);
}

void PDE::step_bc_drainage_contamination()
{
	//  External contamination boundary condition: fixed boundary.
	data->inside[x_len - 1] = i1 * (data->old_inside[x_len - 2] 
		+ data->bag_concentration) + data->old_inside[x_len - 1] 
		* (i2 - i3 * data->old_inside[x_len - 1]) 
		+ f9 * data->get_flow(x_len - 1, r_len - 2);
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
	// Record flow profile as well?
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
	double l = 1; // Length of bladder in contact with bladder
	double w = 5e-3; // Depth to which we consider the bladder volume in contact with the catheter
	double S_c = 2 * 3.14 * param->external_catheter_radius * l; // Surface contact area
	double V_c = 3.14 * (l + w) * (param->external_catheter_radius + w) // Contact volume
		* (param->external_catheter_radius + w) - 3.14 * l 
		* param->external_catheter_radius * param->external_catheter_radius; 
	x_len = data->x_len;
	r_len = data->r_len;
	N = int(time / dt);
	print_step = int(print_interval / dt);
	dx = param->catheter_length / (double(x_len)-1);
	dr = param->catheter_radius / (double(r_len) - 1);
	R_square = std::pow(param->catheter_radius, 2);
	dr_square = dr * dr;
	o1 = param->diffusivity * dt / (dx*dx);
	o2 = 1.0 - 2.0 * o1 + param->growth_rate1 * dt;
	o3 = param->growth_rate1 * dt / param->carrying_capacity1;
	o4 = param->attachment_rate * V_c * dt / (2*3.14*param->external_catheter_radius*dx);
	o5 = o2 - param->detachment_rate * S_c * dt / (2 * 3.14 * param->external_catheter_radius * dx);
	coi = param->diffusivity * dt / (dx * dx);
	b1 = 1.0 + param->growth_rate2 * dt - dt * (param->attachment_rate * V_c + param->urine_rate) / param->sump_volume;
	b2 = param->growth_rate2 * dt / param->carrying_capacity2;
	b3 = param->detachment_rate * S_c * dt / param->sump_volume;
	i1 = param->diffusivity * dt / (dx*dx); // THIS IS EXACTLY THE SAME AS o1
	i2 = 1.0 - 2.0 * i1 + param->growth_rate3*dt;
	i3 = param->growth_rate3 * dt / param->carrying_capacity3;
	f1 = param->diffusivity * dt / (dr_square); 
	f2 = 2 * param->urine_rate * dt / (3.14 * std::pow(param->catheter_radius, 4) * dx); 
	f3 = 1 - 2 * f1; 
	f7b = f2 * R_square; 
	f8b = f3 - f7b; 
	f9 = param->diffusivity * dt/ dr; 
}

void PDE::stability_check()
{
	if (o1 > 0.5) // Stability check for external surface
	{
		std::cout << "D dt / dx^2 = " << o1 << " > 0.5, simulation unstable \n";
		exit(EXIT_FAILURE);
	}
	else if (f1 + f7b > 0.5) // Stability check for inside flow
	{
		std::cout << "( D dt / dr^2  = " << f1 << " ) + (u dt / dx = " << f7b << " ) > 0.5, simulation unstable \n";
		exit(EXIT_FAILURE);
	}
}