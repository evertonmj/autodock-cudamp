/*
 * Wrapper for selection allocation
 * Compiled with Cuda compiler.
 */

#include <cuda.h>
#include <cuda_runtime.h>
#include "cutil.h"
#include <stdio.h>
#include "typedefs.h"
#include "select_alloc_wrapper.h"

#define BLOCK_SIZE 128

/**
 * GPU kernel for alloc
 * @param alloc_data alloc array to be stored into
 */
__global__ void alloc_kernel(float *alloc_data)
{
    int idx = blockIdx.x  * blockDim.x + threadIdx.x;

//    alloc_data[idx] = make_float4(1.0,1.0,1.0,1.0);
    alloc_data[idx] = 1.0;

}

////////////////////////////////////////////////////////////////////////////////
//! Entry point for Cuda function
//! @param num_individuals number of individuals in population
//! @param alloc CPU alloc array
//! @param bratwurst worst value
//! @param energy
//! @param invdiffwa
////////////////////////////////////////////////////////////////////////////////
//extern "C" void
extern "C" void select_alloc_wrapper(unsigned int num_individuals, Real *alloc, double bratwurst, double energy, double indvwa)
{
    const unsigned int mem_size = sizeof(float) * num_individuals;
    cudaError_t kerr;
    int nBlocks;

    // allocate device memory
    float* alloc_data;
    cudaMalloc((void**) &alloc_data, mem_size);

    // setup execution parameters
    nBlocks = num_individuals/BLOCK_SIZE + ((num_individuals%BLOCK_SIZE==0)?0:1);

    // execute the kernel
    alloc_kernel<<< nBlocks, BLOCK_SIZE >>>((float *) alloc_data);
    kerr = cudaGetLastError();

    if (kerr != cudaSuccess)
    { fprintf(stderr, "CUDA ERROR = %s\n", cudaGetErrorString(kerr)); }

    // copy results from device to host
    cudaMemcpy(alloc, alloc_data, mem_size, cudaMemcpyDeviceToHost);

    // cleanup memory
    cudaFree(alloc_data);
}
