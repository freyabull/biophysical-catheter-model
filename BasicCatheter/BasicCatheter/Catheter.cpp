// This structure contains the data at the current timestep, and should be 
// initialized with the appropriate initial conditions for the problem.

#include "pch.h"
#include "Catheter.h"


Catheter::Catheter() :
	bladder(0.0), skin_concentration(0.0),
	bag_concentration(0.0), x_len(0), r_len(0)
{
}

Catheter::Catheter(double t_skin_conc, double t_bag_conc, int t_x_len, 
	               int t_r_len) :
	outside(t_x_len, 0.0), old_outside(t_x_len, 0.0), inside(t_x_len, 0.0), 
	old_inside(t_x_len, 0.0), flow(t_x_len*t_r_len, 0.0), 
	old_flow(t_x_len* t_r_len, 0.0), bladder(0.0), 
	skin_concentration(t_skin_conc), bag_concentration(t_bag_conc), 
	x_len(t_x_len), r_len(t_r_len)
{
}

Catheter::Catheter(std::vector<double> t_outside, double t_bladder,
	               std::vector<double> t_inside, double t_skin_conc, 
	               double t_bag_conc, int t_r_len) :
	outside(t_outside), old_outside(t_outside), bladder(t_bladder), 
	inside(t_inside), old_inside(t_inside), 
	flow(t_outside.size()* t_r_len, t_bladder),
	old_flow(t_outside.size()* t_r_len, t_bladder), 
	skin_concentration(t_skin_conc), bag_concentration(t_bag_conc), 
	x_len(t_outside.size()), r_len(t_r_len)
{
	// Ensure bacterial concentration in flow is zero at surface
	for (int i = 0; i < x_len; ++i) { set_flow(i, 0, 0.0); }
}

Catheter::~Catheter()
{
}

void Catheter::update()
{
	// Swap the old and new data ready for the next timestep.
	std::swap(outside, old_outside);
	std::swap(inside, old_inside);
	std::swap(flow, old_flow);
	// Avoid issues arising from a zero carrying capacity
	if (isnan(bladder)) { bladder = 0.0; } 
}

double Catheter::get_flow(int i, int j)
{
	return old_flow[i * r_len + j];
}

void Catheter::set_flow(int i, int j, double val)
{
	flow[i * r_len + j] = val;
}
