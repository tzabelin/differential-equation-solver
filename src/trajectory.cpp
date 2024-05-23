#include <math.h>
#include "trajectory.h"


void KerrTrajectory::initial_conditions(double r_0, double b_y, double b_z, double theta_init, double* y)
{
	double p[4];

	double theta0 = theta_init;
	double r0 = 2 * r_0;
	double rho2 = pow(r0, 2) + pow(a, 2) * pow(cos(theta0), 2);
	double rho = sqrt(rho2);
	double Delta = pow(r0, 2) - 2.0*r0 + pow(a, 2);
	double by = 2.0*b_y;
	double bz = 2.0*b_z;


	p[0] = 1.0 / rho * ((pow(r0, 2) + pow(a, 2)) / sqrt(Delta) + a * by*sin(theta0) / sqrt(pow(by, 2) + bz * 82 + pow(r0, 2)));
	p[1] = -1 * sqrt(Delta) / rho * r0 / sqrt(pow(by, 2) + pow(bz, 2) + pow(r0, 2));
	p[2] = 1.0 / rho * bz / sqrt(pow(by, 2) + pow(bz, 2) + pow(r0, 2));
	p[3] = 1.0 / rho * (a / sqrt(Delta) + by / (sin(theta0)*sqrt(pow(by, 2) + pow(bz, 2) + pow(r0, 2))));
	double Energy = p[0] * (1.0 - 2 * r0 / rho2) + p[3] * 2.0 * a* r0* pow((sin(theta0)), 2) / rho2;

	y[8] = (p[3] * pow(sin(theta0), 2) * (pow(r0, 2) + pow(a, 2) + 2.0*r0*pow(a, 2) * pow(sin(theta0), 2) / rho2) - p[0] * 2.0*r0*a*pow(sin(theta0), 2) / rho2) / Energy;
	y[7] = pow((rho2*p[2] / Energy), 2) + pow(cos(theta0), 2) * (pow(y[8], 2) / pow(sin(theta0), 2) - pow(a, 2));

	double phi_initial = 0.0;
	y[0] = double(0.0); //set independent variable to 0
	y[1] = 1.0 / r0;
	y[2] = 1.0 + (pow(a, 2) - pow(y[8], 2) - y[7]) / pow(r0, 2) + 2.0*(pow((a - y[8]), 2) + y[7]) / pow(r0, 3) - pow(a, 2) * y[7] / pow(r0, 4);
	if (y[2] < pow(0.1, 8))
	{
		y[2] = 0.0;
	}
	else
	{
		y[2] = sqrt(y[2]);
	}
	y[3] = theta0;
	y[4] = y[7] + pow(cos(theta0), 2) * (pow(a, 2) - pow(y[8], 2) / pow(sin(theta0), 2));

	if (p[2] > 0)
	{
		y[4] = -1 * sqrt(y[4]);
	}

	y[5] = phi_initial;
	y[6] = 0.0;

}

	
