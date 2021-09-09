#include "cuda_runtime.h"
#include "device_launch_parameters.h"

template<typename System, typename precision>
__device__ void RK4(System system, precision *args, const precision dx, const int n, const int args_number)
{
	precision *args_tmp = new precision[args_number];

	// k is different derevatives in different points for calculation
	precision *k1 = new precision[args_number-1];
	precision *k2 = new precision[args_number-1];
	precision *k3 = new precision[args_number-1];
	precision *k4 = new precision[args_number-1];

	system(args, k1);

	args_tmp[0] = args[0] + dx / 2;
	for (int i = 1; i < args_number; i++)
	{
		args_tmp[i] = args[i] + dx / 2 * k1[i-1];
	}
	system(args_tmp, k2);

	args_tmp[0] = args[0] + dx / 2;
	for (int i = 1; i < args_number; i++)
	{
		args_tmp[i] = args[i] + dx/ 2 * (k2[i-1]);
	}
	system(args_tmp, k3);

	args_tmp[0] = args[0] + dx;
	for (int i = 1; i < args_number; i++)
	{
		args_tmp[i] = args[i] + dx * (k3[i-1]);
	}
	system(args_tmp, k4);

	args[n*args_number] = args[0] + dx;
	for (int i = 1; i < args_number; i++)
	{
		args[n*args_number+i] = args[i] + dx / 6 * k1[i-1] + dx / 3 * k2[i-1] + dx / 3 * k3[i-1] + dx / 6 * k4[i-1];
	}

	delete[] args_tmp;
	delete[] k1;
	delete[] k2;
	delete[] k3;
	delete[] k4;
}

