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
	const int n = 55;
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
trajectory.initial_conditions(250, 0.5, 2.3, 3.141592653589793238 / 2, &(args[args_number*4]));
trajectory.initial_conditions(1800, -0.8, -1.5, 3.141592653589793238 / 2, &(args[args_number*5]));
trajectory.initial_conditions(500, 1.2, -2.0, 3.141592653589793238 / 2, &(args[args_number*6]));
trajectory.initial_conditions(700, -1.7, 0.7, 3.141592653589793238 / 2, &(args[args_number*7]));
trajectory.initial_conditions(1200, 0.3, 0.9, 3.141592653589793238 / 2, &(args[args_number*8]));
trajectory.initial_conditions(800, -0.6, -1.8, 3.141592653589793238 / 2, &(args[args_number*9]));
trajectory.initial_conditions(600, 1.5, 2.5, 3.141592653589793238 / 2, &(args[args_number*10]));
trajectory.initial_conditions(1100, -1.1, 1.7, 3.141592653589793238 / 2, &(args[args_number*11]));
trajectory.initial_conditions(350, 0.9, -0.3, 3.141592653589793238 / 2, &(args[args_number*12]));
trajectory.initial_conditions(900, -0.2, 1.4, 3.141592653589793238 / 2, &(args[args_number*13]));
trajectory.initial_conditions(400, 1.8, -1.2, 3.141592653589793238 / 2, &(args[args_number*14]));
trajectory.initial_conditions(1600, -1.5, 0.1, 3.141592653589793238 / 2, &(args[args_number*15]));
trajectory.initial_conditions(1000, 0.7, 1.3, 3.141592653589793238 / 2, &(args[args_number*16]));
trajectory.initial_conditions(1400, -0.4, -2.3, 3.141592653589793238 / 2, &(args[args_number*17]));
trajectory.initial_conditions(1300, 1.6, 0.5, 3.141592653589793238 / 2, &(args[args_number*18]));
trajectory.initial_conditions(1700, -1.4, 1.1, 3.141592653589793238 / 2, &(args[args_number*19]));
trajectory.initial_conditions(1350, 0.1, -1.6, 3.141592653589793238 / 2, &(args[args_number*20]));
trajectory.initial_conditions(950, -0.9, 2.0, 3.141592653589793238 / 2, &(args[args_number*21]));
trajectory.initial_conditions(750, 1.3, -0.5, 3.141592653589793238 / 2, &(args[args_number*22]));
trajectory.initial_conditions(2200, -1.2, 0.3, 3.141592653589793238 / 2, &(args[args_number*23]));
trajectory.initial_conditions(1750, 0.2, -0.7, 3.141592653589793238 / 2, &(args[args_number*24]));
trajectory.initial_conditions(850, -1.6, 1.9, 3.141592653589793238 / 2, &(args[args_number*25]));
trajectory.initial_conditions(1150, 0.8, -1.4, 3.141592653589793238 / 2, &(args[args_number*26]));
trajectory.initial_conditions(1650, -0.5, 1.2, 3.141592653589793238 / 2, &(args[args_number*27]));
trajectory.initial_conditions(1050, 1.4, -0.9, 3.141592653589793238 / 2, &(args[args_number*28]));
trajectory.initial_conditions(1450, -0.3, 2.0, 3.141592653589793238 / 2, &(args[args_number*29]));
trajectory.initial_conditions(1850, 0.6, -1.1, 3.141592653589793238 / 2, &(args[args_number*30]));
trajectory.initial_conditions(950, -1.5, 1.8, 3.141592653589793238 / 2, &(args[args_number*31]));
trajectory.initial_conditions(1350, 0.4, -1.3, 3.141592653589793238 / 2, &(args[args_number*32]));
trajectory.initial_conditions(1750, -0.2, 1.5, 3.141592653589793238 / 2, &(args[args_number*33]));
trajectory.initial_conditions(1950, 0.9, -1.6, 3.141592653589793238 / 2, &(args[args_number*34]));
trajectory.initial_conditions(1250, -1.1, 1.7, 3.141592653589793238 / 2, &(args[args_number*35]));
trajectory.initial_conditions(1650, 0.7, -0.8, 3.141592653589793238 / 2, &(args[args_number*36]));
trajectory.initial_conditions(2050, -0.4, 2.1, 3.141592653589793238 / 2, &(args[args_number*37]));
trajectory.initial_conditions(1450, 1.1, -1.0, 3.141592653589793238 / 2, &(args[args_number*38]));
trajectory.initial_conditions(1850, -0.6, 1.6, 3.141592653589793238 / 2, &(args[args_number*39]));
trajectory.initial_conditions(1550, 0.3, -1.2, 3.141592653589793238 / 2, &(args[args_number*40]));
trajectory.initial_conditions(1950, -0.8, 1.3, 3.141592653589793238 / 2, &(args[args_number*41]));
trajectory.initial_conditions(2150, 1.2, -1.7, 3.141592653589793238 / 2, &(args[args_number*42]));
trajectory.initial_conditions(1050, -0.9, 1.4, 3.141592653589793238 / 2, &(args[args_number*43]));
trajectory.initial_conditions(1550, 0.5, -1.5, 3.141592653589793238 / 2, &(args[args_number*44]));
trajectory.initial_conditions(1950, -1.0, 1.9, 3.141592653589793238 / 2, &(args[args_number*45]));
trajectory.initial_conditions(2250, 1.0, -1.8, 3.141592653589793238 / 2, &(args[args_number*46]));
trajectory.initial_conditions(1150, -0.7, 1.6, 3.141592653589793238 / 2, &(args[args_number*47]));
trajectory.initial_conditions(1750, 0.8, -1.3, 3.141592653589793238 / 2, &(args[args_number*48]));
trajectory.initial_conditions(2050, -0.5, 1.7, 3.141592653589793238 / 2, &(args[args_number*49]));
trajectory.initial_conditions(1350, 1.3, -1.2, 3.141592653589793238 / 2, &(args[args_number*50]));
trajectory.initial_conditions(1950, -0.6, 1.8, 3.141592653589793238 / 2, &(args[args_number*51]));
trajectory.initial_conditions(1850, 0.7, -1.5, 3.141592653589793238 / 2, &(args[args_number*52]));
trajectory.initial_conditions(2150, -0.8, 1.4, 3.141592653589793238 / 2, &(args[args_number*53]));
trajectory.initial_conditions(2250, 0.9, -1.6, 3.141592653589793238 / 2, &(args[args_number*54]));



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
