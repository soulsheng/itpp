
#include "ldpc_bp_decode.cuh"
#include "ldpc_bp_decode_kernel.cuh"

#include <cuda_runtime.h>
#include <thrust/reduce.h>
#include <thrust/device_vector.h>

bool syndrome_check_gpu(int *LLR, int nvar, 
	int* sumX2, int ncheck, 
	int* V, int nmaxX2 ) 
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
	cudaMalloc( (void**)&d_V, ncheck * nmaxX2 * sizeof(int) );
	cudaMemcpy( d_V, V, ncheck * nmaxX2 * sizeof(int), cudaMemcpyHostToDevice );

	dim3 block( 256 );
	dim3 grid( (ncheck + block.x - 1) / block.x );

	syndrome_check_kernel<<< grid, block >>>( d_LLR, d_sumX2, ncheck, d_V, d_synd );

	int sum = thrust::reduce( thrust::device_ptr<int>( d_synd ),
		thrust::device_ptr<int>( d_synd + ncheck ), 
		(int) 0, thrust::plus<int>());

	cudaFree( d_synd );
	cudaFree( d_LLR );
	cudaFree( d_sumX2 );
	cudaFree( d_V );

	return sum == ncheck;   // codeword is valid
}

void updateVariableNode_gpu( int nvar, int ncheck, int nmaxX1, int nmaxX2, 
	int* sumX1, int* mcv, int* mvc, int* iind, int * LLRin, int * LLRout ) 
{

	int* d_sumX1 ;
	cudaMalloc( (void**)&d_sumX1, nvar * sizeof(int) );
	cudaMemcpy( d_sumX1, sumX1, nvar * sizeof(int), cudaMemcpyHostToDevice );
	
	int* d_mcv ;
	cudaMalloc( (void**)&d_mcv, ncheck * nmaxX2 * sizeof(int) );
	cudaMemcpy( d_mcv, mcv, ncheck * nmaxX2 * sizeof(int), cudaMemcpyHostToDevice );
		
	int* d_mvc ;
	cudaMalloc( (void**)&d_mvc, nvar * nmaxX1 * sizeof(int) );
	cudaMemcpy( d_mvc, mvc, nvar * nmaxX1 * sizeof(int), cudaMemcpyHostToDevice );

	int* d_iind ;
	cudaMalloc( (void**)&d_iind, nvar * nmaxX1 * sizeof(int) );
	cudaMemcpy( d_iind, iind, nvar * nmaxX1 * sizeof(int), cudaMemcpyHostToDevice );

	int* d_LLRin ;
	cudaMalloc( (void**)&d_LLRin, nvar * sizeof(int) );
	cudaMemcpy( d_LLRin, LLRin, nvar * sizeof(int), cudaMemcpyHostToDevice );

	int* d_LLRout ;
	cudaMalloc( (void**)&d_LLRout, nvar * sizeof(int) );
	cudaMemcpy( d_LLRout, LLRout, nvar * sizeof(int), cudaMemcpyHostToDevice );

	dim3 block( 256 );
	dim3 grid( (nvar + block.x - 1) / block.x );

	updateVariableNode_kernel<<< grid, block >>>( nvar, d_sumX1, d_mcv, d_mvc, d_iind, d_LLRin, d_LLRout );
	
	cudaMemcpy( LLRout, d_LLRout, nvar * sizeof(int), cudaMemcpyDeviceToHost );
	cudaMemcpy( mvc, d_mvc, nvar * nmaxX1 * sizeof(int), cudaMemcpyDeviceToHost );

	cudaFree( d_sumX1 );
	cudaFree( d_mcv );
	cudaFree( d_mvc );
	cudaFree( d_iind );
	cudaFree( d_LLRin );
	cudaFree( d_LLRout );

}

void updateCheckNode_gpu( int nvar, int ncheck, int nmaxX1, int nmaxX2, 
	int* sumX2, int* mcv, int* mvc, int* jind, 
	short int Dint1, short int Dint2, short int Dint3, int* logexp_table,
	int* jj, int* m, int* ml, int* mr, int max_cnd, int QLLR_MAX )
{
	
	int* d_sumX2 ;
	cudaMalloc( (void**)&d_sumX2, ncheck * sizeof(int) );
	cudaMemcpy( d_sumX2, sumX2, ncheck * sizeof(int), cudaMemcpyHostToDevice );

	int* d_mcv ;
	cudaMalloc( (void**)&d_mcv, ncheck * nmaxX2 * sizeof(int) );
	cudaMemcpy( d_mcv, mcv, ncheck * nmaxX2 * sizeof(int), cudaMemcpyHostToDevice );
		
	int* d_mvc ;
	cudaMalloc( (void**)&d_mvc, nvar * nmaxX1 * sizeof(int) );
	cudaMemcpy( d_mvc, mvc, nvar * nmaxX1 * sizeof(int), cudaMemcpyHostToDevice );

	int* d_jind ;
	cudaMalloc( (void**)&d_jind, nvar * nmaxX2 * sizeof(int) );
	cudaMemcpy( d_jind, jind, ncheck * nmaxX2 * sizeof(int), cudaMemcpyHostToDevice );
	
	int* d_logexp_table ;
	cudaMalloc( (void**)&d_logexp_table, Dint2 * sizeof(int) );
	cudaMemcpy( d_logexp_table, logexp_table, Dint2 * sizeof(int), cudaMemcpyHostToDevice );
	
	int* d_jj ;
	cudaMalloc( (void**)&d_jj, max_cnd * sizeof(int) );
	cudaMemset( d_jj, 0, max_cnd * sizeof(int) );
	
	int* d_m ;
	cudaMalloc( (void**)&d_m, max_cnd * sizeof(int) );
	cudaMemset( d_m, 0, max_cnd * sizeof(int) );
	
	int* d_ml ;
	cudaMalloc( (void**)&d_ml, max_cnd * sizeof(int) );
	cudaMemset( d_ml, 0, max_cnd * sizeof(int) );
	
	int* d_mr ;
	cudaMalloc( (void**)&d_mr, max_cnd * sizeof(int) );
	cudaMemset( d_mr, 0, max_cnd * sizeof(int) );

	dim3 block( 256 );
	dim3 grid( (ncheck + block.x - 1) / block.x );

	updateCheckNode_kernel<<< grid, block >>>(ncheck, 
		d_sumX2, d_mcv, d_mvc, d_jind, Dint1, Dint2, Dint3, d_logexp_table,
		d_jj, d_m, d_ml, d_mr, QLLR_MAX );

	cudaMemcpy( mcv, d_mcv, ncheck * nmaxX2 * sizeof(int), cudaMemcpyDeviceToHost );
	cudaMemcpy( mvc, d_mvc, nvar * nmaxX1 * sizeof(int), cudaMemcpyDeviceToHost );

	cudaFree( d_sumX2 );
	cudaFree( d_mcv );
	cudaFree( d_mvc );
	cudaFree( d_jind );
	cudaFree( d_logexp_table );
	cudaFree( d_jj );	cudaFree( d_m );	cudaFree( d_ml );	cudaFree( d_mr );
	
}