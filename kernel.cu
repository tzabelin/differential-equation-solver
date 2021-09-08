#include <iostream>
#include <stdio.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "RK.cu"

struct  Simple
{
	__device__ void operator()(const double *x, const double *y, double *k, const int &n)
	{
		*k = sin(*x);
	}
};

int main()
{
	const int n = 2;
	const int size = 500;

	double *x = new double[n*size];
	x[0] = double(1.0);
	x[1] = double(4.0);
	double *y = new double[n*size];
	y[0] = double(1.0);
	y[1] = double(16.0);
	double h(1.0);

	double*x_dev;
	double*y_dev;
	double*h_dev;

	cudaSetDevice(0);

	if ((cudaMalloc((void**)&x_dev, sizeof(double)*n*size) ||
		cudaMalloc((void**)&y_dev, sizeof(double)*n*size) ||
		cudaMalloc((void**)&h_dev, sizeof(h))) != cudaSuccess)
	{
		std::cerr << "CUDA malloc error";
		exit(-1);
	}

	if ((cudaMemcpy(x_dev, x, sizeof(double)*n*size, cudaMemcpyHostToDevice) ||
		cudaMemcpy(y_dev, y, sizeof(double)*n*size, cudaMemcpyHostToDevice) ||
		cudaMemcpy(h_dev, &h, sizeof(h), cudaMemcpyHostToDevice)) != cudaSuccess)
	{
		std::cerr << "CUDA mempcy error";
		exit(-1);
	}

	Simple s;
	do_step << <1, 2 >> > (s, x_dev, y_dev, h_dev, n, size);
	cudaDeviceSynchronize();

	if ((
		cudaMemcpy(x, x_dev, sizeof(double)*n*size, cudaMemcpyDeviceToHost) ||
		cudaMemcpy(y, y_dev, sizeof(double)*n*size, cudaMemcpyDeviceToHost)) != cudaSuccess)
	{
		std::cerr << "CUDA mempcy device to host error";
		exit(-1);
	}

	cudaFree(x_dev);
	cudaFree(y_dev);
	cudaFree(h_dev);

	for (int u = 0; u < n*size; u += 2)
	{
		std::cout << x[u];
		std::cout << "    ";
		std::cout << y[u];
		std::cout << "    ";
		std::cout << y[u + 1];
		std::cout << "    ";
		std::cout << std::endl;
	}
}