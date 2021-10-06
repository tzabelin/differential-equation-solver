#pragma once
#include "cuda_runtime.h"
#include <math.h>


struct  KerrTrajectory
{

	const double a = 0.9982;

	/*operator () defines system of differential equations*/
	__device__ void KerrTrajectory::operator()(const double *args, double *k)
	{
		k[0] = args[2];
		k[1] = -2 * a * a * args[7] * args[1] * args[1] * args[1] + 3 * ((a - args[8])*(a - args[8]) + args[7]) * args[1] * args[1] + (a * a - args[8] * args[8] - args[7]) * args[1];
		k[2] = args[4];
		k[3] = cos(args[3]) * (args[8] * args[8] / pow(sin(args[3]), 3) - a * a * sin(args[3]));
		k[4] = -1.0 * (a - args[8] / pow(sin(args[3]), 2)) + a * (1 / (args[1] * args[1]) + a * a - args[8] * a) / (1 / (args[1] * args[1]) - 2 / args[1] + a * a);
		k[5] = -1.0 * a * (a * sin(args[3]) * sin(args[3]) - args[8]) + (1 / (args[1] * args[1]) + a * a) * (1 / (args[1] * args[1]) + a * a - args[8] * a) / (1 / (args[1] * args[1]) - 2 / args[1] + a * a);
		k[6] = 0;
		k[7] = 0;
	}

	/*sets initial values to variables in system according to given distance, deflection on two axes and angle*/
	void initial_conditions(double r_0, double b_y, double b_z, double theta_init, double* y);
	
};