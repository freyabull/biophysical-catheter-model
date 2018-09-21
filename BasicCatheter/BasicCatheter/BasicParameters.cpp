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
	double t_droplet_size, double t_stickiness, double t_sump_volume,
	double t_catheter_length) :
	diffusivity(t_diffusivity), growth_rate1(t_growth_rate1),
	carrying_capacity1(t_carrying_capacity1), growth_rate2(t_growth_rate2),
	carrying_capacity2(t_carrying_capacity2), growth_rate3(t_growth_rate3),
	carrying_capacity3(t_carrying_capacity3), urine_rate(t_urine_rate),
	droplet_size(t_droplet_size), stickiness(t_stickiness),
	sump_volume(t_sump_volume), catheter_length(t_catheter_length)
{
}

BasicParameters::~BasicParameters()
{
}
