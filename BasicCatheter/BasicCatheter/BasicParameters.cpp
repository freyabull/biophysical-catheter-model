// This structure contains the parameters relevant to the physics of the problem.

#include "pch.h"
#include "BasicParameters.h"


BasicParameters::BasicParameters()
{
}

BasicParameters::BasicParameters(double t_diffusivity, double t_growth_rate1,
	double t_carrying_capacity1, double t_growth_rate2,
	double t_carrying_capacity2, double t_growth_rate3,
	double t_carrying_capacity3, double t_urine_rate,
	double t_catheter_radius, double t_stickiness, double t_sump_volume,
	double t_catheter_length, double t_coupling_o, double t_coupling_i) :
	diffusivity(t_diffusivity), growth_rate1(t_growth_rate1),
	carrying_capacity1(t_carrying_capacity1), growth_rate2(t_growth_rate2),
	carrying_capacity2(t_carrying_capacity2), growth_rate3(t_growth_rate3),
	carrying_capacity3(t_carrying_capacity3), urine_rate(t_urine_rate),
	catheter_radius(t_catheter_radius), stickiness(t_stickiness),
	sump_volume(t_sump_volume), catheter_length(t_catheter_length),
	coupling_o(t_coupling_o), coupling_i(t_coupling_i)
{
}

BasicParameters::~BasicParameters()
{
}
