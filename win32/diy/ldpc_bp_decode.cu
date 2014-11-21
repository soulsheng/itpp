
#include "ldpc_bp_decode.cuh"
#include "ldpc_bp_decode_kernel.cuh"

#include <cuda_runtime.h>
#include <thrust/reduce.h>
#include <thrust/device_vector.h>

bool syndrome_check_gpu(int *LLR, int nvar, 
	int* sumX2, int ncheck, 
	int* V, int maxX2 ) 
{
	// Please note the IT++ convention that a sure zero corresponds to
	// LLR=+infinity
	int* d_synd ;
	cudaMalloc( (void**)&d_synd, ncheck * sizeof(int) );
	cudaMemset( d_synd, 0, ncheck * sizeof(int) );

	int* d_LLR ;
	cudaMalloc( (void**)&d_LLR, nvar * sizeof(int) );
	cudaMemcpy( d_LLR, LLR, nvar * sizeof(int), cudaMemcpyHostToDevice );

	int* d_sumX2 ;
	cudaMalloc( (void**)&d_sumX2, ncheck * sizeof(int) );
	cudaMemcpy( d_sumX2, sumX2, ncheck * sizeof(int), cudaMemcpyHostToDevice );

	int* d_V ;
	cudaMalloc( (void**)&d_V, ncheck * maxX2 * sizeof(int) );
	cudaMemcpy( d_V, V, ncheck * maxX2 * sizeof(int), cudaMemcpyHostToDevice );

	dim3 block( 256 );
	dim3 grid( (ncheck + block.x - 1) / block.x );

	syndrome_check_kernel<<< grid, block >>>( d_LLR, d_sumX2, ncheck, d_V, d_synd );

	int bValid = thrust::reduce( thrust::device_ptr<int>( d_synd ),
		thrust::device_ptr<int>( d_synd + ncheck ), 
		(int) 0, thrust::multiplies<int>());

	cudaFree( d_synd );
	cudaFree( d_LLR );
	cudaFree( d_sumX2 );
	cudaFree( d_V );

	return bValid;   // codeword is valid
}
