
#include "driver-updateVar.cuh"
#include "ldpc_bp_decode_kernel.cuh"
#include "driverUtility.h"

#include <cuda_runtime.h>


bool driverUpdataVar::launch()
{
	dim3 block( SIZE_BLOCK );
	dim3 grid( (nvar + block.x - 1) / block.x );

	updateVariableNode_kernel<<< grid, block >>>( nvar, ncheck, 
		d_sumX1, d_mcv, d_iind, d_input, 
		d_output, d_mvc );

	cudaError_t	status = cudaGetLastError();
	return status == cudaSuccess ;
}

bool driverUpdataVar::verify()
{
	cudaMemcpy( mvc, d_mvc, nvar * nmaxX1 * sizeof(int), cudaMemcpyHostToDevice );
	cudaMemcpy( output, d_output, nvar * sizeof(char), cudaMemcpyDeviceToHost );

	// mvc
	int i = 0;
	for ( ; i < nvar * nmaxX1; i++ )
	{
		if ( ref_mvc[i] != mvc[i] )
			break;
	}

	if ( i < nvar * nmaxX1 )
		return false;

	// output
	int j = 0;
	for ( ; j < nvar; j++ )
	{
		if ( ref_output[j] != output[j] )
			break;
	}

	if ( j < nvar )
		return false;

	return true;
}

driverUpdataVar::driverUpdataVar()
	: nvar( VAR_SIZE_CODE )
	, ncheck( CHECK_SIZE_CODE )
	, nmaxX1( MAX_VAR_NODE )
	, nmaxX2( MAX_CHECK_NODE )
{
	sumX1 = (int*)malloc(nvar * sizeof(int));
	sumX2 = (int*)malloc(ncheck * sizeof(int));
	d_iind = (int*)malloc(nvar * nmaxX1 * sizeof(int));
	d_jind = (int*)malloc(ncheck * nmaxX2 * sizeof(int));
	mvc = (int*)malloc(nvar * nmaxX1 * sizeof(int));
	mcv = (int*)malloc(ncheck * nmaxX2 * sizeof(int));
	input = (int*)malloc(nvar * sizeof(int));
	output = (char*)malloc(nvar * sizeof(char));

	ref_mvc = (int*)malloc(nvar * nmaxX1 * sizeof(int));
	ref_output = (char*)malloc(nvar * sizeof(char));

	readArray( sumX1, nvar, "../data/sumX1.txt" );
	readArray( sumX2, ncheck, "../data/sumX2.txt" );

	readArray( iind, nvar * nmaxX1, "../data/iind.txt" );
	readArray( jind, ncheck * nmaxX2, "../data/jind.txt" );

	readArray( ref_output, nvar, "../data/output.txt" );
	readArray( ref_mvc, nvar * nmaxX1, "../data/mvc.txt" );		

	cudaMalloc( (void**)&d_sumX1, nvar * sizeof(int) );		// const 64 K
	cudaMemcpy( d_sumX1, sumX1, nvar * sizeof(int), cudaMemcpyHostToDevice );

	cudaMalloc( (void**)&d_sumX2, ncheck * sizeof(int) );	// const 32 K
	cudaMemcpy( d_sumX2, sumX2, ncheck * sizeof(int), cudaMemcpyHostToDevice );

	cudaMalloc( (void**)&d_iind, nvar * nmaxX1 * sizeof(int) );		// const 1.2 M
	cudaMemcpy( d_iind, iind, nvar * nmaxX1 * sizeof(int), cudaMemcpyHostToDevice );

	cudaMalloc( (void**)&d_jind, ncheck * nmaxX2 * sizeof(int) );	// const 300 K
	cudaMemcpy( d_jind, jind, ncheck * nmaxX2 * sizeof(int), cudaMemcpyHostToDevice );

	cudaMalloc( (void**)&d_mcv, ncheck * nmaxX2 * sizeof(int) );
	cudaMemcpy( d_mcv, mcv, ncheck * nmaxX2 * sizeof(int), cudaMemcpyHostToDevice );

	cudaMalloc( (void**)&d_mvc, nvar * nmaxX1 * sizeof(int) );
	cudaMemcpy( d_mvc, mvc, nvar * nmaxX1 * sizeof(int), cudaMemcpyHostToDevice );

	cudaMalloc( (void**)&d_input, nvar * sizeof(int) );
	cudaMalloc( (void**)&d_output, nvar * sizeof(char) );

}

driverUpdataVar::~driverUpdataVar()
{
	// host
	free(sumX1);	free(sumX2);
	free(iind);		free(jind);
	free(mvc);		free(mcv);
	free(input);	free(output);

	free(ref_mvc);	free(ref_output);

	// device
	cudaFree( d_sumX1 );	cudaFree( d_sumX2 );
	cudaFree( d_iind );		cudaFree( d_jind );
	cudaFree( d_mvc );		cudaFree( d_mcv );
	cudaFree( d_input );	cudaFree( d_output );
}