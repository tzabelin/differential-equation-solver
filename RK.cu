#include "cuda_runtime.h"
#include "device_launch_parameters.h"

template<typename System, typename precision = double>
__device__ void RK4(System system, precision *x, precision *y, const precision *dx, const int n, const int step)
{
	precision dx2 = (*dx) / 2;
	precision dx3 = (*dx) / 3;
	precision dx6 = (*dx) / 6;

	precision *y_tmp = new precision(0.0);
	precision *x_tmp = new precision(0.0);

	precision *k1 = new precision(0.0);
	precision *k2 = new precision(0.0);
	precision *k3 = new precision(0.0);
	precision *k4 = new precision(0.0);

	system(x, y_tmp, k1, n);

	*y_tmp = *y + dx2 * (*k1);
	*x_tmp = *x + dx2;
	system(x_tmp, y_tmp, k2, n);

	*y_tmp = *y + dx2 * (*k2);
	*x_tmp = *x + dx2;
	system(x_tmp, y_tmp, k3, n);

	*y_tmp = *y + (*dx) * (*k3);
	*x_tmp = *x + (*dx);
	system(x_tmp, y_tmp, k4, n);

	y[n] = precision(0.0);
	x[n] = precision(0.0);
	y[n] = (*y) + dx6 * (*k1) + dx3 * (*k2) + dx3 * (*k3) + dx6 * (*k4);
	x[n] = (*x) + (*dx);
}
template<typename System, typename precision = double>
__global__ void do_step(System system, precision *x, precision *y, const precision *dx, const int n, const int size)
{
	for (int i = 0; i < size; i++)
	{
		RK4(system, &(x[i*n + threadIdx.x]), &(y[i*n + threadIdx.x]), dx, n, i);
	}
}
