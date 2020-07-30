// This structure contains the parameters relevant to the physics of the problem.

#include "pch.h"
#include "BasicParameters.h"


BasicParameters::BasicParameters() :
	diffusivity(0.0), growth_rate1(0.0), carrying_capacity1(0.0), 
	growth_rate2(0.0), carrying_capacity2(0.0), growth_rate3(0.0),
	carrying_capacity3(0.0), urine_rate(0.0), catheter_radius(0.0), 
	external_catheter_radius(0.0), stickiness(0.0), sump_volume(0.0), 
	catheter_length(0.0), attachment_rate(0.0), detachment_rate(0.0)
{
}

BasicParameters::BasicParameters(double t_diffusivity, double gr1, double cc1, 
	double gr2, double cc2, double gr3, double cc3, double t_urine_rate,
	double in_radius, double ex_radius, double t_stickiness, double volume,
	double length, double attach, double detach) :
	diffusivity(t_diffusivity), growth_rate1(gr1), carrying_capacity1(cc1), 
	growth_rate2(gr2), carrying_capacity2(cc2), growth_rate3(gr3), 
	carrying_capacity3(cc3), urine_rate(t_urine_rate), 
	catheter_radius(in_radius), external_catheter_radius(ex_radius), 
	stickiness(t_stickiness), sump_volume(volume), catheter_length(length),
	attachment_rate(attach), detachment_rate(detach)
{
}

BasicParameters::~BasicParameters()
{
}
