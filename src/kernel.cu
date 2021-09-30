#include <iostream>
#include <stdio.h>
#include<math.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "RK.cu"
#include "trajectory.h"


template<typename System, typename precision = double>
__global__ void do_step(System system, precision *args, const precision dx, const int n, const int size, const int args_number)
{
	for (int i = 0; i < size*n*args_number; i+=n*args_number)
	{
		Runge_Kutta::Order_2(system, &(args[i + threadIdx.x*args_number]), dx, n, args_number);
	}
}


int main()
{
	//number of systems of equations
	const int n = 2;
	//number of points to be processed
	const int size = 1200;
	//number of variable in one system, including independent
	const int args_number = 9;

	double *args = new double[n*size*args_number];

	//Starting conditions
	KerrTrajectory trajectory;
	trajectory.initial_conditions(100, -1.3, 2.587146695207, 3.141592653589793238 / 2, args);
	trajectory.initial_conditions(150, 1.0, 1.0, 3.141592653589793238 / 2, &(args[9]));

	const double h(0.001);

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

	do_step << <1, n >> > (trajectory, args_dev, h, n, size, args_number);
	cudaDeviceSynchronize();

	if (cudaMemcpy(args, args_dev, sizeof(double)*n*size*args_number, cudaMemcpyDeviceToHost) != cudaSuccess)
	{
		std::cerr << "CUDA mempcy device to host error";
		exit(-1);
	}

	cudaFree(args_dev);

	double R = 0, x=0, y=0, z=0;
	for (int i = 0; i < n*size*args_number; i=i+args_number*n)
	{
		R = 1 / (args[i+1] * 2);
		x =R * cos(args[i+3]) * sin(args[i+5]);
		y =R * sin(args[i+3]) * sin(args[i + 5]);
		z =R * cos(args[i + 5]);
		std::cout << x;
		std::cout << "    ";
		std::cout << y;
		std::cout << "    ";
		std::cout << z;
		std::cout << "    ";

		R = 1 / (args[i + 10] * 2);
		x = R * cos(args[i + 12]) * sin(args[i + 5]);
		y = R * sin(args[i + 12]) * sin(args[i + 5]);
		z = R * cos(args[i + 14]);
		std::cout << x;
		std::cout << "    ";
		std::cout << y;
		std::cout << "    ";
		std::cout << z;
		std::cout << "    ";

		std::cout << std::endl;
	}
	delete[] args;
}