
#include "driver-updateVar.cuh"

#include <cuda_runtime.h>

#include <iostream>
using namespace std;


#define		MAX_CHECK_NODE		10
#define		MAX_VAR_NODE		19
#define		VAR_SIZE_CODE		16200
#define		CHECK_SIZE_CODE		8073
#define		SIZE_BLOCK			256

__global__ 
void updateVariableNodeOpti_kernel( const int nvar, const int ncheck, const int* sumX1, const int* mcv, const int* iind, const int * LLRin, 
	char * LLRout, int* mvc ) // not used, just for testing performance bound
{	//	mcv const(input)-> mvc (output)
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	
	if( i>= nvar )
		return;
	

	int mvc_temp = LLRin[i];

	int m[MAX_VAR_NODE];

	for (int jp = 0; jp < sumX1[i]; jp++)
		m[jp] = mcv[ iind[i + jp*nvar] ];


	for (int jp = 0; jp < sumX1[i]; jp++)
		mvc_temp += m[jp];
	

	LLRout[i] = mvc_temp<0;
	
	for (int jp = 0; jp < sumX1[i]; jp++)
			mvc[i + jp*nvar] = mvc_temp - m[jp];

}

bool driverUpdataVar::launch()
{
	dim3 block( SIZE_BLOCK );
	dim3 grid( (nvar + block.x - 1) / block.x );

	updateVariableNodeOpti_kernel<<< grid, block >>>( nvar, ncheck, 
		d_sumX1, d_mcv, d_iind, d_input, 
		d_output, d_mvc );

	cudaError_t	status = cudaGetLastError();
	return status == cudaSuccess ;
}

bool driverUpdataVar::verify()
{
	cudaMemcpy( mvc, d_mvc, nvar * nmaxX1 * sizeof(int), cudaMemcpyDeviceToHost );
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

template <typename T>
void 	readArray(T* pArray, int nSize, char* strFileName)
{
	FILE* fp = NULL;
	fp = fopen( strFileName, "rb" );
	if(fp == NULL)
	{
		printf("failed to open: %s!\n", strFileName);
	}
	fread( pArray, sizeof(T), nSize, fp);
	fclose(fp);
}

driverUpdataVar::driverUpdataVar()
	: nvar( VAR_SIZE_CODE )
	, ncheck( CHECK_SIZE_CODE )
	, nmaxX1( MAX_VAR_NODE )
	, nmaxX2( MAX_CHECK_NODE )
{
	sumX1 = (int*)malloc(nvar * sizeof(int));
	sumX2 = (int*)malloc(ncheck * sizeof(int));
	iind = (int*)malloc(nvar * nmaxX1 * sizeof(int));
	jind = (int*)malloc(ncheck * nmaxX2 * sizeof(int));
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

	readArray( input, nvar, "../data/input.txt" );

	readArray( mcv, ncheck * nmaxX2, "../data/mcv.txt" );	

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

	cudaMalloc( (void**)&d_input, nvar * sizeof(int) );
	cudaMemcpy( d_input, input, nvar * sizeof(int), cudaMemcpyHostToDevice );
	
	cudaMalloc( (void**)&d_output, nvar * sizeof(char) );
	cudaMemset( d_output, 0, nvar * sizeof(char) );

	cudaMalloc( (void**)&d_mvc, nvar * nmaxX1 * sizeof(int) );
	cudaMemset( d_mvc, 0, nvar * nmaxX1 * sizeof(int) );

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