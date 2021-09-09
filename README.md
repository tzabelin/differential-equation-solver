# differential-equation-solver
 Simple solver for system of differential equations in parallel. It produces array of dots that can be visualized using different plotters such as gnuplot.

# Prerequisites

 1)Git
 
 2)NVIDIA CUDA 10 or higher

# Installation

 run *nvcc kernel.cu RK.cu* with appropriate *--arch* and *-gencode* in the src folder
