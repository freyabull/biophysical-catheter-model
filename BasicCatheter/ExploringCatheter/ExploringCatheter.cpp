// An alternate main function that simulates multiple catheters to enable an 
// exploration of parameter space.

#include "../BasicCatheter/Catheter.h"
#include "../BasicCatheter/BasicParameters.h"
#include "../BasicCatheter/PDE.h"
#include <sstream>

int main() {
	// File root for saving
	std::string file_root = "results/urine_rates";

	// Parameter being explored
	// dilution rate: physically relevant range is from 25/6 mm^3 s^-1 to 2500/6 mm^3 s^-1
	double urine_rates[6] = { 5,10,20,30,40,50 };

	// Define variables (bar exploration variable)
	double diffusivity = 1e-2; // diffusivity of bacteria on catheter, in mm^2/s
	double growth_rate1 = std::log(2.0) / 3600.0; // Growth rate of bacteria on outside of catheter, per s
	double carrying_capacity1 = 1e6; // Carrying capacity of bacteria on outside of catheter, per mm^2
	double growth_rate2 = std::log(2.0) / 1800.0; // Growth rate of bacteria in bladder, per s
	double carrying_capacity2 = 1e8; // Carrying capacity of bacteria in bladder, per mm^3 
	double growth_rate3 = std::log(2.0) / 3600.0; // Growth rate of bacteria on inside of catheter, per s
	double carrying_capacity3 = 1e6; // Carrying capacity of bacteria on inside of catheter, per mm^2
	//double urine_rate = 1.0e3 / 60.0; // Rate of generation of urine by kidneys, in mm^3/s 
	double catheter_radius = 1.0; // Radius of catheter in mm
	double stickiness = 1.0; // Probability of a bacterium sticking to surface if it comes into contact
	double sump_volume = 5.0e4; // Volume of residual urine in bladder, in mm^3
	double catheter_length = 40;//160.0; // Length of catheter, in mm
	double skin_concentration = 1e2; // Concentration of bacteria on skin (ie boundary condition)
	double bag_concentration = -0.1; // Concentration of bacteria in drainage bag (ie boundary condition)
	int initial_condition = 1; // Type of initial conditions
	int x_len = 11;// 101; // Number of x points
	int r_len = 11; // Number of r points
	double initial_skin_conc = 1e2; // Concentration of bacteria on the skin at time of catheter insertion
	double lambda = 1.0 / catheter_length; // Scale factor for exponential distributions
	int simulation_length = 86400 * 2; // Timeframe of simulation, in s
	double dt = 0.1; // Time step
	int print_interval = 3600; // Time interval at which to output data (s)
	double catheter_external_radius = catheter_radius;// +1.0; // External catheter radius in mm
	double attachment_rate = 4 * 3.14 * diffusivity * 1e-3; // Rate at which a bacterium in contact sticks (s^-1) Smoluchowski - 4 pi D sigma
	double detachment_rate = 0.05; // Rate at which a bacterium in contact detaches (s^-1)

	// Loop over exploration variable
	for (int i=0; i < sizeof(urine_rates)/sizeof(double); i++) {
		// Set exploration parameter value
		double urine_rate = urine_rates[i];
		// Open file
		std::ofstream results_file; // File for output to be written to 
		// Kind of silly work around to cast i to str type
		std::ostringstream ss;
		ss << i;
		std::string file_name = file_root + "_" + ss.str() + ".csv";
		std::cout << file_name << std::endl;
		results_file.open(file_name, std::ios::trunc);
		if (!results_file.is_open()) {
			std::cout << "Unable to open file " << file_name << std::endl;
			break;
		}
		// Set up catheter
		BasicParameters myParam = BasicParameters(diffusivity, growth_rate1, // Parameters needed to solve the PDE problem
			carrying_capacity1, growth_rate2, carrying_capacity2, growth_rate3,
			carrying_capacity3, urine_rate, catheter_radius, catheter_external_radius, stickiness,
			sump_volume, catheter_length, attachment_rate, detachment_rate);
		Catheter myCatheter = Catheter(skin_concentration, bag_concentration, x_len, r_len); // Current state of catheter
		PDE myPDE = PDE(&myParam, &myCatheter, simulation_length, dt, print_interval);
		// Run simulation
		results_file << "simulation length,time step,print interval,diffusivity,outside growth rate,"
			<< "outside carrying capacity,bladder growth rate,bladder carrying capacity," <<
			"inside growth rate,inside carrying capacity,urine rate,catheter radius,external catheter radius,stickiness,sump volume,"
			<< "catheter length,initial condition,num of x steps,skin concentration," <<
			"drainage bag concentration,attachment rate,detachment rate" << "\n";
		results_file << simulation_length << "," << dt << "," << print_interval << "," << diffusivity <<
			"," << growth_rate1 << "," << carrying_capacity1 << "," << growth_rate2 << "," <<
			carrying_capacity2 << "," << growth_rate3 << "," << carrying_capacity3 << "," << urine_rate << "," <<
			catheter_radius << "," << catheter_external_radius << "," << stickiness << "," << sump_volume << "," << catheter_length <<
			"," << initial_condition << "," << x_len << "," << skin_concentration << "," <<
			bag_concentration << "," << attachment_rate << "," << detachment_rate << "\n\n";
		myPDE.solve(results_file);
		// Close file
		results_file.close();
	}
		
}