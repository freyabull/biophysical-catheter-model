#pragma once
#include <vector>
#include <functional>

class Integrate1D // Implementation of 1D numerical integration methods
{
public: 
	// Structure containing result and error
	struct results { double sol, err; };
	// Trapezoidal integrator: takes as arguments function, bounds, N, returns result struct
	static results trapezoidal(std::function<double (double)> func, double a, double b, int N=100);
	// Discrete trapezoidal integrator: takes as arguments array of y values, x bounds, returns result struct
	static results discrete_trapezoidal(std::vector<double> f, double a, double b);
};


