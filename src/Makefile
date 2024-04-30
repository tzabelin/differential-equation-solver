# Makefile for CUDA project
NVCC = nvcc
CC = g++
CUDA_PATH = /usr/local/cuda
NVCCFLAGS = -I$(CUDA_PATH)/include
CXXFLAGS = -I$(CUDA_PATH)/include
LDFLAGS = -L$(CUDA_PATH)/lib64 -lcudart

# Specify the target file and the install directory
output: trajectory.o kernel.o RK.o
	$(CC) -o output trajectory.o kernel.o RK.o $(LDFLAGS)

kernel.o: kernel.cu
	$(NVCC) $(NVCCFLAGS) -c kernel.cu

RK.o: RK.cu
	$(NVCC) $(NVCCFLAGS) -c RK.cu

trajectory.o: trajectory.cpp trajectory.h
	$(CC) $(CXXFLAGS) -c trajectory.cpp

clean:
	rm -f *.o output
