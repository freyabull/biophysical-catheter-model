// An alternate main function that simulates multiple catheters to enable an 
// exploration of parameter space.

#include "../BasicCatheter/Catheter.h"
#include "../BasicCatheter/BasicParameters.h"
#include "../BasicCatheter/PDE.h"
#include <sstream>

int main() {
	// File root for saving
	std::string file_root = "results/";

	// Parameter being explored
	double catheter_lengths[9] = { 40,60,80,100,120,140,160,180,200 };

	// Define variables (bar exploration variable)
	double diffusivity = 1e-4; // diffusivity of bacteria on catheter, in mm^2/s
	double surface_diffusivity = 1e-4; // diffusivity of bacteria on catheter, in mm^2/s
	double growth_rate1 = std::log(2.0) / 3600.0; // Growth rate of bacteria on outside of catheter, per s
	double carrying_capacity1 = 1e7; // Carrying capacity of bacteria on outside of catheter, per mm^2
	double growth_rate2 = std::log(2.0) / 1800.0; // Growth rate of bacteria in bladder, per s
	double carrying_capacity2 = 1e6; // Carrying capacity of bacteria in bladder, per mm^3 
	double growth_rate3 = std::log(2.0) / 3600.0; // Growth rate of bacteria on inside of catheter, per s
	double carrying_capacity3 = 1e7; // Carrying capacity of bacteria on inside of catheter, per mm^2
	double urine_rate = 1.0e3 / 60.0; // Rate of generation of urine by kidneys, in mm^3/s 
	double catheter_radius = 1.0; // Radius of catheter in mm
	double stickiness = 1.0; // Probability of a bacterium sticking to surface if it comes into contact
	double sump_volume = 5.0e4; // Volume of residual urine in bladder, in mm^3
	double catheter_length = 40.0; // Length of catheter, in mm
	double skin_concentration = 1e2; // Concentration of bacteria on skin (ie boundary condition)
	double bag_concentration = -0.1; // Concentration of bacteria in drainage bag (ie boundary condition)
	int initial_condition = 1; // Type of initial conditions
	int x_len = catheter_length * 3;// 170;// 250; // Number of x points
	int print_num_steps = 101; // Number of x steps at which to output data
	int r_len = 11; // Number of r points
	//double initial_skin_conc = 1e2; // Concentration of bacteria on the skin at time of catheter insertion
	int simulation_length = 86400 * 14;// 200; // Timeframe of simulation, in s
	double dt = 60;// 80; // Time step
	int print_interval = 600; // Time interval at which to output data (s)
	double catheter_external_radius = catheter_radius;// +1.0; // External catheter radius in mm
	double attachment_rate = 4 * 3.14 * diffusivity * 1e-3; // Rate at which a bacterium in contact sticks (s^-1) Smoluchowski - 4 pi D sigma
	double detachment_rate = growth_rate1; // Rate at which a bacterium in contact detaches (s^-1)
	double viscosity = 0.83;

	// Loop over exploration variable
		for (int j = 0; j < sizeof(catheter_lengths) / sizeof(double); j++) {
			// Set exploration parameter value
			double catheter_length = catheter_lengths[j];
			// Open file
			std::ofstream results_file; // File for output to be written to 
			// Kind of silly work around to cast i to str type
			std::ostringstream ss;
			int k = j;
			if (k < 10) { ss << 0; }
			if (k < 100) { ss << 0; }
			if (k < 1000) { ss << 0; }
			ss << k;
			std::string file_name = file_root + ss.str() + ".csv";
			std::cout << file_name << std::endl;
			results_file.open(file_name, std::ios::trunc);
			if (!results_file.is_open()) {
				std::cout << "Unable to open file " << file_name << std::endl;
				break;
			}
			// Set up catheter
			BasicParameters myParam = BasicParameters(diffusivity, surface_diffusivity, growth_rate1, // Parameters needed to solve the PDE problem
				carrying_capacity1, growth_rate2, carrying_capacity2, growth_rate3,
				carrying_capacity3, urine_rate, catheter_radius, catheter_external_radius, stickiness,
				sump_volume, catheter_length, viscosity, attachment_rate, detachment_rate);
			Catheter myCatheter = Catheter(skin_concentration, bag_concentration, x_len, r_len); // Current state of catheter
			PDE myPDE = PDE(&myParam, &myCatheter, simulation_length, dt, print_interval, print_num_steps);
			// Run simulation
			results_file << "simulation length,time step,print interval,print num steps,diffusivity,surface diffusivity,outside growth rate,"
				<< "outside carrying capacity,bladder growth rate,bladder carrying capacity," <<
				"inside growth rate,inside carrying capacity,urine rate,catheter radius,external catheter radius,stickiness,sump volume,"
				<< "catheter length,initial condition,num of x steps,skin concentration," <<
				"drainage bag concentration,viscosity,attachment rate,detachment rate" << "\n";
			results_file << simulation_length << "," << dt << "," << print_interval << "," << print_num_steps << "," << diffusivity << "," << surface_diffusivity <<
				"," << growth_rate1 << "," << carrying_capacity1 << "," << growth_rate2 << "," <<
				carrying_capacity2 << "," << growth_rate3 << "," << carrying_capacity3 << "," << urine_rate << "," <<
				catheter_radius << "," << catheter_external_radius << "," << stickiness << "," << sump_volume << "," << catheter_length <<
				"," << initial_condition << "," << x_len << "," << skin_concentration << "," <<
				bag_concentration << "," << viscosity << "," << attachment_rate << "," << detachment_rate << "\n\n";
			myPDE.solve_light(results_file);
			// Close file
			results_file.close();
		}
		
}