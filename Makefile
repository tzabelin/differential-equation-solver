# Makefile for CUDA project
NVCC = nvcc
CC = g++
CUDA_PATH = /usr/local/cuda
NVCCFLAGS = -I$(CUDA_PATH)/include
CXXFLAGS = -I$(CUDA_PATH)/include
LDFLAGS = -L$(CUDA_PATH)/lib64 -lcudart
SRC_DIR = src
BUILD_DIR = build

# Specify the target file and the install directory
output: $(BUILD_DIR)/trajectory.o $(BUILD_DIR)/kernel.o $(BUILD_DIR)/RK.o
	$(CC) -o output $(BUILD_DIR)/trajectory.o $(BUILD_DIR)/kernel.o $(BUILD_DIR)/RK.o $(LDFLAGS)

$(BUILD_DIR)/kernel.o: $(SRC_DIR)/kernel.cu
	$(NVCC) $(NVCCFLAGS) -c $(SRC_DIR)/kernel.cu

$(BUILD_DIR)/RK.o: $(SRC_DIR)/RK.cu
	$(NVCC) $(NVCCFLAGS) -c $(SRC_DIR)/RK.cu

$(BUILD_DIR)/trajectory.o: $(SRC_DIR)/trajectory.cpp $(SRC_DIR)/trajectory.h
	$(CC) $(CXXFLAGS) -c $(SRC_DIR)/trajectory.cpp

clean:
	rm -f $(BUILD_DIR)/*.o output
