// Implementation of numerical integration methods

#include "pch.h"
#include "Integrate.h"

Integrate1D::results Integrate1D::trapezoidal(std::function<double(double)> func, double a,
	double b, int N)
{
	double dx = (b - a) / N; // Grid spacing
	double result = 0.5 * ( func(a) + func(b) );
	for (int i = 1; i < N; i++) {
		result += func(a + i * dx);
	}
	result = result * dx;
	double error = -(b - a) * (b - a) / (12 * N * N * dx) * (func(b) 
		- func(b - dx) - func(a + dx) + func(a));
	return results{ result, error };
}

Integrate1D::results Integrate1D::discrete_trapezoidal(std::vector<double> f, 
	double a, double b) 
{
	int N = f.size(); // Number of grid points
	double dx = (b - a) / (N-1); // Grid spacing
	double result = 0.5 * (f[0] + f[N - 1]);
	for (int i = 1; i < N - 1; i++) {
		result += f[i];
	}
	result = result * dx;
	double error = -(b - a) * (b - a) / (12 * N * N * dx) * (f[N-1] - f[N-2]
		- f[1] + f[0]);
	return results{ result, error };
}