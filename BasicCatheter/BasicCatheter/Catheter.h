#pragma once
#include <vector>
#include <utility>
class Catheter // Data for the current time-step, should be initialised with any initial and/or boundary conditions
{
public:
	std::vector<double> outside; // Concentrations of bacteria on outside of catheter
	std::vector<double> old_outside; // Concentrations of bacteria on outside of catheter for previous time-step
	std::vector<double> inside; // Concentrations of bacteria on inside of catheter
	std::vector<double> old_inside; // Concentrations of bacteria on inside of catheter for previous time-step
	double bladder; // Concentration of bacteria within bladder
	double skin_concentration; // Concentration of bacteria on the skin, where a negative value = "switched off"

	Catheter();
	Catheter(double skin_concentration, int x_len);
	Catheter(std::vector<double> outside, double bladder, 
		     std::vector<double> inside, double skin_concentration = -1.0);
	~Catheter();

	void update();
};

