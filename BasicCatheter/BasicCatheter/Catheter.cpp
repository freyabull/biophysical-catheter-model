// This structure contains the data at the current timestep, and should be 
// initialized with the appropriate initial conditions for the problem.

#include "pch.h"
#include "Catheter.h"


Catheter::Catheter()
{
}

Catheter::Catheter(double t_skin_conc, double t_bag_conc, int x_len) :
	outside(x_len, 0.0), old_outside(x_len, 0.0), inside(x_len, 0.0), 
	old_inside(x_len, 0.0), bladder(0.0), skin_concentration(t_skin_conc),
	bag_concentration(t_bag_conc), x_len(x_len)
{
}

Catheter::Catheter(std::vector<double> t_outside, double t_bladder,
	               std::vector<double> t_inside, double t_skin_conc, 
	               double t_bag_conc) :
	outside(t_outside), old_outside(t_outside), bladder(t_bladder), 
	inside(t_inside), old_inside(t_inside), skin_concentration(t_skin_conc),
	bag_concentration(t_bag_conc), x_len(t_outside.size())
{
}

Catheter::~Catheter()
{
}

void Catheter::update()
{
	// Swap the old and new data ready for the next timestep.
	std::swap(outside, old_outside);
	std::swap(inside, old_inside);
	if (isnan(bladder)) { bladder = 0.0; } // Avoid issues arising from a zero carrying capacity
}
