
#include "ldpc_bp_decode.cuh"
#include "ldpc_bp_decode_kernel.cuh"
//#include "driverUtility.h"

#include <cuda_runtime.h>
#include <thrust/reduce.h>
#include <thrust/device_vector.h>

#if USE_TEXTURE_ADDRESS
	cudaArray* arr_mcv;
	cudaArray* arr_mvc;
	cudaChannelFormatDesc channelDesc;
#endif

bool ldpc_gpu::syndrome_check_gpu() 
{
	dim3 block( SIZE_BLOCK );
	dim3 grid( (ncheck + block.x - 1) / block.x );

	syndrome_check_kernel<<< grid, block >>>( d_LLRout, d_sumX2, ncheck, d_V, d_synd );

	int h_synd=0;
	cudaMemcpy( &h_synd, d_synd, sizeof(int), cudaMemcpyDeviceToHost );

	return h_synd == 0;   // codeword is valid
}

void ldpc_gpu::updateVariableNode_gpu() 
{
	dim3 block( SIZE_BLOCK );
	dim3 grid( (nvar + block.x - 1) / block.x );

	updateVariableNode_kernel<<< grid, block >>>( nvar, ncheck, d_sumX1, d_mcv, d_iind, d_LLRin, d_LLRout, d_mvc );
}

void ldpc_gpu::updateCheckNode_gpu()
{
	dim3 block( SIZE_BLOCK );
	dim3 grid( (ncheck + block.x - 1) / block.x );

	updateCheckNode_kernel<<< grid, block >>>(ncheck, nvar, 
		d_sumX2, d_mvc, d_jind, d_logexp_table, Dint1, Dint2, Dint3,
		QLLR_MAX, d_mcv );	// Shared not faster
}

void ldpc_gpu::initializeMVC_gpu( )
{
	dim3 block( 256 );
	dim3 grid( (nvar + block.x - 1) / block.x );

	initializeMVC_kernel<<< grid, block >>>( nvar, d_sumX1, d_LLRin, d_mvc );
}

int ldpc_gpu::bp_decode(int *LLRin, int *LLRout,
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

#if USE_TEXTURE_ADDRESS
    // update the array to the texture
    cudaMemcpyToArray(arr_mcv, 0, 0, d_mcv, ncheck * nmaxX2 * sizeof(int), cudaMemcpyDeviceToDevice);
#endif

    // step 2: variable to check nodes
	updateVariableNode_gpu();

#if USE_TEXTURE_ADDRESS
    // update the array to the texture
    cudaMemcpyToArray(arr_mvc, 0, 0, d_mvc, nvar * nmaxX1 * sizeof(int), cudaMemcpyDeviceToDevice);
#endif

#if	USE_TABLE_CODE
	updateConstantMemoryLLRByte( d_LLRout );
#endif

	if (psc && syndrome_check_gpu()) {
	  is_valid_codeword = true;
      break;
    }
  }
  while (iter < max_iters);

  cudaMemcpy( LLRout, d_LLRout, nvar * sizeof(int), cudaMemcpyDeviceToHost );


  return (is_valid_codeword ? iter : -iter);
}

int ldpc_gpu::bp_decode_once(int *LLRin, char *LLRout,
	bool psc /*= true*/,			//!< check syndrom after each iteration
	int max_iters /*= 50*/ )		//!< Maximum number of iterations
{
	cudaMemcpy( d_LLRin, LLRin, nvar * sizeof(int), cudaMemcpyHostToDevice );

 	dim3 block( SIZE_BLOCK );
	dim3 grid( (nvar + block.x - 1) / block.x );

	// initial step
	initializeMVC_kernel<<< grid, block >>>( nvar, d_sumX1, d_LLRin, d_mvc );

#if WRITE_FILE_FOR_DRIVER
	static bool bRunOnce1 = false;
	if( !bRunOnce1 ){
		cudaMemcpy( h_mvc, d_mvc, nvar * nmaxX1 * sizeof(int), cudaMemcpyDeviceToHost );
		writeArray( h_mvc, nvar * nmaxX1, "../data/mvcInit.txt" );		
		bRunOnce1 = true;
	}
#endif

	int not_valid_codeword = true;
	int iter = 1;
	for( ; iter < max_iters && not_valid_codeword; iter ++ )
	{
		// --------- Step 1: check to variable nodes ----------
		updateCheckNodeOpti_kernel<<< grid, block >>>(ncheck, nvar, 
			d_sumX2, d_mvc, d_jind, d_logexp_table, Dint1, Dint2, Dint3,QLLR_MAX, 
			d_mcv );	// Shared not faster

				
#if WRITE_FILE_FOR_DRIVER
		static bool bRunOnce1 = false;
		if( iter == 1 && !bRunOnce1 ){
			cudaMemcpy( h_mcv, d_mcv, ncheck * nmaxX2 * sizeof(int), cudaMemcpyDeviceToHost );

			writeArray( h_mcv, ncheck * nmaxX2, "../data/mcv.txt" );

			bRunOnce1 = true;
		}
#endif

		// --------- Step 2: variable to check nodes ----------
#if USE_BLOCK_2D
	
	dim3 block( SIZE_BLOCK_2D_X, MAX_VAR_NODE );
	dim3 grid;
	grid.x = (nvar * MAX_VAR_NODE + SIZE_BLOCK_2D_X * MAX_VAR_NODE - 1) 
				/ (SIZE_BLOCK_2D_X * MAX_VAR_NODE) ;

	updateVariableNodeOpti2D_kernel<<< grid, block >>>( nvar, ncheck, 
		d_sumX1, d_mcv, d_iind, d_LLRin, 
		d_LLRout, d_mvc );
#else
		updateVariableNodeOpti_kernel<<< grid, block >>>( nvar, ncheck, 
			d_sumX1, d_mcv, d_iind, d_LLRin, 
			d_LLRout, d_mvc );
#endif

		// --------- Step 3: check syndrome ��żУ�� ----------
#if 0
		syndrome_check_kernel<<< grid, block >>>( d_LLRout, d_sumX2, ncheck, d_V, 
			d_synd );

		cudaMemcpy( &not_valid_codeword, d_synd, sizeof(int), cudaMemcpyDeviceToHost );
#else
		cudaMemcpy( LLRout, d_LLRout, nvar * sizeof(char), cudaMemcpyDeviceToHost );
		
#if WRITE_FILE_FOR_DRIVER
		static bool bRunOnce2 = false;
		if( iter == 1 && !bRunOnce2 ){
			cudaMemcpy( h_mvc, d_mvc, nvar * nmaxX1 * sizeof(int), cudaMemcpyDeviceToHost );

			writeArray( LLRout, nvar, "../data/output.txt" );
			writeArray( h_mvc, nvar * nmaxX1, "../data/mvc.txt" );		

			bRunOnce2 = true;
		}
#endif

		if (psc && check_parity_cpu(LLRout)) {
			 not_valid_codeword = false;
			break;
		}
#endif
	}
  
  //cudaMemcpy( LLRout, d_LLRout, nvar * sizeof(char), cudaMemcpyDeviceToHost );


  return (!not_valid_codeword ? iter : -iter);
}

bool ldpc_gpu::check_parity_cpu(char *LLR) 
{
	// Please note the IT++ convention that a sure zero corresponds to
	// LLR=+infinity
	int i, j, synd, vi;

	for (j = 0; j < ncheck; j++) {
		synd = 0;
		int vind = j; // tracks j+i*ncheck
		for (i = 0; i < h_sumX2[j]; i++) {
			vi = h_V[vind];
			if (LLR[vi]) {
				synd++;
			}
			vind += ncheck;
		}
		if ((synd&1) == 1) {
			return false;  // codeword is invalid
		}
	}
	return true;   // codeword is valid
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
	
	this->h_V = V;
	this->h_sumX2 = sumX2;

	//max_cnd = 200;
	QLLR_MAX = (std::numeric_limits<int>::max() >> 4);

	cudaMalloc( (void**)&d_LLRin, nvar * sizeof(int) );
	cudaMalloc( (void**)&d_LLRout, nvar * sizeof(char) );
	cudaMemset( d_LLRout, 1, nvar * sizeof(char) );

	cudaMalloc( (void**)&d_synd, 1 * sizeof(int) );
	cudaMemset( d_synd, 0, 1 * sizeof(int) );
	
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
	cudaMemset( d_mcv, 0, ncheck * nmaxX2 * sizeof(int) );
		
	cudaMalloc( (void**)&d_mvc, nvar * nmaxX1 * sizeof(int) );
	cudaMemset( d_mvc, 0, nvar * nmaxX1 * sizeof(int) );

	cudaMalloc( (void**)&d_logexp_table, Dint2 * sizeof(int) );		// const 1.2 K
	cudaMemcpy( d_logexp_table, logexp_table, Dint2 * sizeof(int), cudaMemcpyHostToDevice );

	initConstantMemoryLogExp(logexp_table);

#if USE_TEXTURE_ADDRESS
	// cuda texture ------------------------------------------------------------------------------------------
	channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindSigned);
    cudaError_t err = cudaMallocArray(&arr_mcv, &channelDesc, ncheck, nmaxX2);
    cudaMemcpyToArray(arr_mcv, 0, 0, d_mcv, ncheck * nmaxX2 * sizeof(int), cudaMemcpyDeviceToDevice);

	texMCV.addressMode[0] = cudaAddressModeClamp;
	texMCV.addressMode[1] = cudaAddressModeClamp;
    texMCV.filterMode = cudaFilterModePoint;
    texMCV.normalized = false;

	cudaBindTextureToArray(texMCV, arr_mcv, channelDesc);

	cudaMallocArray(&arr_mvc, &channelDesc, nvar, nmaxX1);
    cudaMemcpyToArray(arr_mvc, 0, 0, d_mvc, nvar * nmaxX1 * sizeof(int), cudaMemcpyDeviceToDevice);
	cudaBindTextureToArray(texMVC, arr_mvc, channelDesc);

#endif

	h_mvc = (int*)malloc(nvar * nmaxX1 * sizeof(int));
	h_mcv = (int*)malloc(ncheck * nmaxX2 * sizeof(int));

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
	
	cudaFree( d_logexp_table );	

	free( h_mvc );	free( h_mcv );

	return true;
}

ldpc_gpu::~ldpc_gpu()
{
	release();
}