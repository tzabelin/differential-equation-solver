#include <iostream>
#include <stdio.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "RK.cu"

struct  Simple
{
	__device__ void operator()(const double *args, double *k)
	{
		k[0] = sin(args[0]);
		k[1] = cos(args[1]);
	}
};
template<typename System, typename precision = double>
__global__ void do_step(System system, precision *args, const precision dx, const int n, const int size, const int args_number)
{
	for (int i = 0; i < size*n*args_number; i+=n*args_number)
	{
		RK4(system, &(args[i + threadIdx.x*args_number]), dx, n, args_number);
	}
}
int main()
{
	//number of systems of equations
	const int n = 2;
	//number of points to be processed
	const int size = 500;
	//number of variable in one system
	const int args_number = 3;

	//Starting conditions
	double *args = new double[n*size*args_number];
	args[0] = double(1.0);
	args[1] = double(1.0);
	args[2] = double(1.0);

	args[3] = double(4.0);
	args[4] = double(9.0);
	args[5] = double(16.0);

	double h(1.0);

	double *args_dev;

	cudaSetDevice(0);

	if (cudaMalloc((void**)&args_dev, sizeof(double)*n*size*args_number) != cudaSuccess)
	{
		std::cerr << "CUDA malloc error";
		exit(-1);
	}

	if (cudaMemcpy(args_dev, args, sizeof(double)*n*size*args_number, cudaMemcpyHostToDevice)  != cudaSuccess)
	{
		std::cerr << "CUDA mempcy error";
		exit(-1);
	}

	Simple s;
	do_step << <1, n >> > (s, args_dev, h, n, size, args_number);
	cudaDeviceSynchronize();

	if (cudaMemcpy(args, args_dev, sizeof(double)*n*size*args_number, cudaMemcpyDeviceToHost) != cudaSuccess)
	{
		std::cerr << "CUDA mempcy device to host error";
		exit(-1);
	}

	cudaFree(args_dev);

	for (int i = 0; i < n*size*args_number; i ++)
	{
		std::cout << args[i];
		std::cout << "    ";
		if ((i+1) % (n*args_number) == 0)
		{
			std::cout << std::endl;
		}
	}
	delete[] args;
}