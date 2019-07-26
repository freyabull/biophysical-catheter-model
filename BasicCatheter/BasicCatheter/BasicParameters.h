#pragma once
class BasicParameters // Parameters defining the physics of the problem
{
public:
	double diffusivity; // Diffusivity of bacteria
	double growth_rate1; // Growth rate on outside of catheter
	double carrying_capacity1; // Carrying capacity in bladder
	double growth_rate2; // Growth rate in bladder
	double carrying_capacity2; // Carrying capacity on inside of catheter
	double growth_rate3; // Growth rate on outside of catheter
	double carrying_capacity3; // Carrying capacity on inside of catheter
	double urine_rate; // Rate of production of urine
	double catheter_radius; // Internal radius of catheter
	double stickiness; // Probability of a bacterium sticking to wall if it comes in contact
	double sump_volume; // Volume of residual urine
	double catheter_length; // Length of the catheter
	double coupling_o = 1.0; // Coupling strength between outside and bladder
	double coupling_i = 1.0; // Coupling strength between bladder and inside
	double coupling_oi = 1.0; // Coupling strength between outside and inside

	BasicParameters();
	BasicParameters(double diffusivity, double growth_rate1,
		double carrying_capacity1, double growth_rate2,
		double carrying_capacity2, double growth_rate3,
		double carrying_capacity3, double urine_rate,
		double catheter_radius, double stickiness, double sump_volume,
		double catheter_length, double coupling_o, double coupling_i,
		double coupling_oi);
	~BasicParameters();
};

