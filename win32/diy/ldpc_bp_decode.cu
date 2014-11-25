
#include "ldpc_bp_decode.cuh"
#include "ldpc_bp_decode_kernel.cuh"

#include <cuda_runtime.h>
#include <thrust/reduce.h>
#include <thrust/device_vector.h>

bool ldpc_gpu::syndrome_check_gpu() 
{
	dim3 block( 256 );
	dim3 grid( (ncheck + block.x - 1) / block.x );

	syndrome_check_kernel<<< grid, block >>>( d_LLRout, d_sumX2, ncheck, d_V, d_synd );

	int sum = thrust::reduce( thrust::device_ptr<int>( d_synd ),
		thrust::device_ptr<int>( d_synd + ncheck ), 
		(int) 0, thrust::plus<int>());

	return sum == ncheck;   // codeword is valid
}

void ldpc_gpu::updateVariableNode_gpu() 
{
	dim3 block( 256 );
	dim3 grid( (nvar + block.x - 1) / block.x );

	updateVariableNode_kernel<<< grid, block >>>( nvar, d_sumX1, d_mcv, d_mvc, d_iind, d_LLRin, d_LLRout );
}

void ldpc_gpu::updateCheckNode_gpu()
{
	dim3 block( 256 );
	dim3 grid( (ncheck + block.x - 1) / block.x );

	updateCheckNode_kernel<<< grid, block >>>(ncheck, 
		d_sumX2, d_mcv, d_mvc, d_jind, Dint1, Dint2, Dint3,
		d_jj, d_m, d_ml, d_mr, max_cnd, QLLR_MAX );
}

void ldpc_gpu::initializeMVC_gpu( )
{
	dim3 block( 256 );
	dim3 grid( (nvar + block.x - 1) / block.x );

	initializeMVC_kernel<<< grid, block >>>( nvar, d_sumX1, d_mvc, d_LLRin );
}

int ldpc_gpu::bp_decode(int *LLRin, int *LLRout,
	int* sumX1,	int* mvc, 
	bool psc /*= true*/,			//!< check syndrom after each iteration
	int max_iters /*= 50*/ )		//!< Maximum number of iterations
{
	cudaMemcpy( d_LLRin, LLRin, nvar * sizeof(int), cudaMemcpyHostToDevice );

  // initial step
	initializeMVC_gpu();

  bool is_valid_codeword = false;
  int iter = 0;
  do {
    iter++;
    //if (nvar >= 100000) { it_info_no_endl_debug("."); }
    // --------- Step 1: check to variable nodes ----------
	updateCheckNode_gpu();

    // step 2: variable to check nodes
	updateVariableNode_gpu();

	if (psc && syndrome_check_gpu()) {
	  is_valid_codeword = true;
      break;
    }
  }
  while (iter < max_iters);

  cudaMemcpy( LLRout, d_LLRout, nvar * sizeof(int), cudaMemcpyDeviceToHost );


  return (is_valid_codeword ? iter : -iter);
}

bool ldpc_gpu::initialize( int nvar, int ncheck,
	int nmaxX1, int nmaxX2,
	int* sumX1, int* sumX2, int* iind, int* jind, int* V, 	// Parity check matrix parameterization
	int* mvc, int* mcv,	// temporary storage for decoder (memory allocated when codec defined)
	short int Dint1, short int Dint2, short int Dint3,
	int* logexp_table		//! The lookup tables for the decoder
	)
{
	this->nvar = nvar;		this->ncheck = ncheck;
	this->nmaxX1 = nmaxX1;	this->nmaxX2 = nmaxX2; // max(sumX1) max(sumX2)
	this->Dint1 = Dint1;	this->Dint2 = Dint2;	this->Dint3 = Dint3;	//! Decoder (lookup-table) parameters
	
	max_cnd = 200;
	QLLR_MAX = (std::numeric_limits<int>::max() >> 4);

	cudaMalloc( (void**)&d_LLRin, nvar * sizeof(int) );
	cudaMalloc( (void**)&d_LLRout, nvar * sizeof(int) );
	cudaMemset( d_LLRout, 0, nvar * sizeof(int) );

	cudaMalloc( (void**)&d_synd, ncheck * sizeof(int) );
	cudaMemset( d_synd, 0, ncheck * sizeof(int) );
	
	cudaMalloc( (void**)&d_sumX1, nvar * sizeof(int) );		// const 64 K
	cudaMemcpy( d_sumX1, sumX1, nvar * sizeof(int), cudaMemcpyHostToDevice );

	cudaMalloc( (void**)&d_sumX2, ncheck * sizeof(int) );	// const 32 K
	cudaMemcpy( d_sumX2, sumX2, ncheck * sizeof(int), cudaMemcpyHostToDevice );

	cudaMalloc( (void**)&d_iind, nvar * nmaxX1 * sizeof(int) );		// const 1.2 M
	cudaMemcpy( d_iind, iind, nvar * nmaxX1 * sizeof(int), cudaMemcpyHostToDevice );
	
	cudaMalloc( (void**)&d_jind, ncheck * nmaxX2 * sizeof(int) );	// const 300 K
	cudaMemcpy( d_jind, jind, ncheck * nmaxX2 * sizeof(int), cudaMemcpyHostToDevice );

	cudaMalloc( (void**)&d_V, ncheck * nmaxX2 * sizeof(int) );		// const 300 K
	cudaMemcpy( d_V, V, ncheck * nmaxX2 * sizeof(int), cudaMemcpyHostToDevice );
	
	cudaMalloc( (void**)&d_mcv, ncheck * nmaxX2 * sizeof(int) );
	cudaMemcpy( d_mcv, mcv, ncheck * nmaxX2 * sizeof(int), cudaMemcpyHostToDevice );
		
	cudaMalloc( (void**)&d_mvc, nvar * nmaxX1 * sizeof(int) );
	cudaMemcpy( d_mvc, mvc, nvar * nmaxX1 * sizeof(int), cudaMemcpyHostToDevice );

	//cudaMalloc( (void**)&d_logexp_table, Dint2 * sizeof(int) );		// const 1.2 K
	//cudaMemcpy( d_logexp_table, logexp_table, Dint2 * sizeof(int), cudaMemcpyHostToDevice );

	initConstantMemory(logexp_table);

	cudaMalloc( (void**)&d_jj, ncheck * max_cnd * sizeof(int) );
	cudaMemset( d_jj, 0, ncheck * max_cnd * sizeof(int) );
	
	cudaMalloc( (void**)&d_m, ncheck * max_cnd * sizeof(int) );
	cudaMemset( d_m, 0, ncheck * max_cnd * sizeof(int) );
	
	cudaMalloc( (void**)&d_ml, ncheck * max_cnd * sizeof(int) );
	cudaMemset( d_ml, 0, ncheck * max_cnd * sizeof(int) );
	
	cudaMalloc( (void**)&d_mr, ncheck * max_cnd * sizeof(int) );
	cudaMemset( d_mr, 0, ncheck * max_cnd * sizeof(int) );

	return true;
}


bool ldpc_gpu::release()
{
	cudaFree( d_LLRin );	cudaFree( d_LLRout );
	
	cudaFree( d_synd );

	cudaFree( d_sumX1 );	cudaFree( d_sumX2 );
	
	cudaFree( d_iind );		cudaFree( d_jind );
	cudaFree( d_V );

	cudaFree( d_mcv );		cudaFree( d_mvc );
	
	//cudaFree( d_logexp_table );	

	cudaFree( d_jj );	cudaFree( d_m );	cudaFree( d_ml );	cudaFree( d_mr );

	return true;
}

ldpc_gpu::~ldpc_gpu()
{
	release();
}