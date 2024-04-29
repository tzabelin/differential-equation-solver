#include <iostream>
#include <stdio.h>
#include<math.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "RK.cu"
#include "trajectory.h"
#include <chrono>


template<typename System, typename precision = double>
__global__ void do_step_RK(System system, precision *args, const precision dx, const int n, const int size, const int args_number)
{
	for (int i = 0; i < size*n*args_number; i+=n*args_number)
	{
		Runge_Kutta::Order_4(system, &(args[i + threadIdx.x*args_number]), dx, n, args_number);
	}
}
template<typename System, typename precision = double>
__global__ void do_step(System system, precision* args, precision* k_prev, const precision dx, const int n, const int size, const int args_number) {

	//Adams_Bashforth::RK_8(system, &(args[ threadIdx.x * args_number]), &(k_prev[threadIdx.x * (args_number - 1)]), dx, n, args_number);
	for (int i = 0; i < size; i += 1)
	{
		
		int base_idx = i * n * args_number + threadIdx.x * args_number;
		int k_base_idx = i * n * (args_number-1) + threadIdx.x * (args_number - 1);

		Adams_Bashforth::Adams_Bashforth_2<System,precision>(system, &args[base_idx], &k_prev[k_base_idx], dx, n, args_number);
	}
}

int main()
{
	//number of systems of equations
	const int n = 4;
	//number of points to be processed
	const int size = 1200;//00
	//number of variable in one system, including independent
	const int args_number = 9;

	double *args = new double[n*size*args_number];
	double *k_prev = new double[n * size * (args_number-1)];// previous derivative for Adams_Bashforth_2


	//Starting conditions
	KerrTrajectory trajectory;
	trajectory.initial_conditions(100, -1.3, 2.587146695207, 3.141592653589793238 / 2, args);
	trajectory.initial_conditions(150, 1.0, 1.0, 3.141592653589793238 / 2, &(args[args_number]));
	trajectory.initial_conditions(2000, -1.0, 1.0, 3.141592653589793238 / 2, &(args[args_number*2]));
	trajectory.initial_conditions(300, 1.0, -1.0, 3.141592653589793238 / 2, &(args[args_number*3]));



	const double h(0.001);

	double *args_dev;
	double* k_prev_dev; // previous derivative for Adams_Bashforth_2
	trajectory.initial_derivatives(args, k_prev);
	trajectory.initial_derivatives(&(args[args_number]), &(k_prev[(args_number - 1)]));

	cudaSetDevice(0);

	if (cudaMalloc((void**)&args_dev, sizeof(double) * n * size * args_number) != cudaSuccess ||
		cudaMalloc((void**)&k_prev_dev, sizeof(double) * n * size * args_number) != cudaSuccess) {
		std::cerr << "CUDA malloc error";
		exit(-1);
	}

	if ( (cudaMemcpy(args_dev, args, sizeof(double)*n*size*args_number, cudaMemcpyHostToDevice)  != cudaSuccess) 
		|| (cudaMemcpy(k_prev_dev, k_prev, sizeof(double) * n * size * (args_number - 1), cudaMemcpyHostToDevice) != cudaSuccess))
	{
		std::cerr << "CUDA mempcy error";
		exit(-1);
	}

	do_step_RK << <1, n >> > (trajectory, args_dev, h, n, size, args_number); //RK
	//do_step << <1, n >> > (trajectory, args_dev, k_prev_dev, h, n, size, args_number); //Adams_Bashforth_2
	cudaDeviceSynchronize();


	if (cudaMemcpy(args, args_dev, sizeof(double)*n*size*args_number, cudaMemcpyDeviceToHost) != cudaSuccess)
	{
		std::cerr << "CUDA mempcy device to host error";
		exit(-1);
	}

	cudaFree(args_dev);

	double R = 0, x=0, y=0, z=0;
	bool *is_feasible = new bool[args_number]; // check that photons dont go into infinity
	for(int i=0; i<args_number; i++)
	{
		is_feasible[i] = true;
	}
	for (int i = 0; i < n * size * args_number; i = i + args_number * n)
	{
		for (int j = 0; j < n; j++) {
			if (is_feasible[j] = true)
			{
				R = 1 / (args[i + 1 + j * args_number] * 2);
				if (R > 5000)
					is_feasible[j] = false;
				x = R * cos(args[i + 3 + j * args_number]) * sin(args[i + 5 + j * args_number]);
				y = R * sin(args[i + 3 + j * args_number]) * sin(args[i + 5 + j * args_number]);
				z = R * cos(args[i + 5 + j * args_number]);
				std::cout << x;
				std::cout << "    ";
				std::cout << y;
				std::cout << "    ";
				std::cout << z;
				std::cout << "    ";
			}
			else
			{
				std::cout << 0;
				std::cout << "    ";
				std::cout << 0;
				std::cout << "    ";
				std::cout << 0;
				std::cout << "    ";
			}
		}
		std::cout << std::endl;
	}
	delete[] args;
}