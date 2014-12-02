
#pragma once

#define		TABLE_SIZE_DINT2	300
#define		MAX_CHECK_NODE		10
#define		TABLE_SIZE_CODE		16200
#define		USE_TABLE_CODE		0
#define		USE_TEXTURE_ADDRESS	0
#define		USE_SHARED_MLR		0
#define		USE_SHARED_MVC		0
#define		SIZE_BLOCK			256

#if USE_TEXTURE_ADDRESS
texture<int, 2, cudaReadModeElementType> texMCV;
texture<int, 2, cudaReadModeElementType> texMVC;
#endif

__device__ __constant__ int const_logexp_table[TABLE_SIZE_DINT2];

#if		USE_TABLE_CODE
__device__ __constant__ char const_llr_byte[TABLE_SIZE_CODE];	// char

void updateConstantMemoryLLRByte(char *bLLR)
{
	cudaMemcpyToSymbol( const_llr_byte, bLLR, TABLE_SIZE_CODE * sizeof(char), 0, cudaMemcpyDeviceToDevice );
}
#endif

void initConstantMemoryLogExp(int *logexp_table)
{
	cudaMemcpyToSymbol( const_logexp_table, logexp_table, TABLE_SIZE_DINT2 * sizeof(int), 0, cudaMemcpyHostToDevice );
}

__global__ 
void syndrome_check_kernel(const int *d_LLR,
	const int* d_sumX2, const int ncheck, 
	const int* d_V,
	int* d_synd) 
{
	int j = blockIdx.x * blockDim.x + threadIdx.x;

	if( j>= ncheck )
		return;

	if( j == 0 )
		*d_synd = 0;

	__syncthreads();

	int i, vi;
	int synd = 0;
	int vind = j; // tracks j+i*ncheck
	for (i = 0; i < d_sumX2[j]; i++) {
		vi = d_V[vind];
#if	USE_TABLE_CODE
		if ( const_llr_byte[vi] ) {
#else
		if (d_LLR[vi] < 0) {
#endif
			synd++;
		}
		vind += ncheck;
	}

	if( synd&1 )
		atomicAdd( d_synd, 1 );	// synd is even ?
}

__global__ 
void updateVariableNode_kernel( const int nvar, const int ncheck, const int* sumX1, const int* mcv, const int* iind, const int * LLRin, 
	int * LLRout, int* mvc, char* bLLR ) 
{	//	mcv const(input)-> mvc (output)
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	
	if( i>= nvar )
		return;

	int mvc_temp = LLRin[i];
	for (int jp = 0; jp < sumX1[i]; jp++) {
		int index = iind[i + jp*nvar];
#if USE_TEXTURE_ADDRESS
		mvc_temp +=  tex2D(texMCV, index%ncheck, index/ncheck);
#else
		mvc_temp +=  mcv[index];
#endif
	}
	LLRout[i] = mvc_temp;
	int index_iind = i;  // tracks i+j*nvar
	for (int j = 0; j < sumX1[i]; j++) {
		mvc[index_iind] = mvc_temp - mcv[iind[index_iind]];
		index_iind += nvar;
	}
	
	bLLR[i] = LLRout[i]<0;
}

__device__
int  logexp_device(const int x,
	const short int Dint1, const short int Dint2, const short int Dint3	//! Decoder (lookup-table) parameters
	)		
{
	int ind = x >> Dint3;
	if (ind >= Dint2) // outside table
		return 0;

	// Without interpolation
	return const_logexp_table[ind];
}

__device__
int Boxplus(const int a, const int b,
	const short int Dint1, const short int Dint2, const short int Dint3,	//! Decoder (lookup-table) parameters
	const int QLLR_MAX )		//! The lookup tables for the decoder
{
	int a_abs = (a > 0 ? a : -a);
	int b_abs = (b > 0 ? b : -b);
	int minabs = (a_abs > b_abs ? b_abs : a_abs);
	int term1 = (a > 0 ? (b > 0 ? minabs : -minabs)
		: (b > 0 ? -minabs : minabs));

	if (Dint2 == 0) {  // logmax approximation - avoid looking into empty table
		// Don't abort when overflowing, just saturate the QLLR
		if (term1 > QLLR_MAX) {
			return QLLR_MAX;
		}
		if (term1 < -QLLR_MAX) {
			return -QLLR_MAX;
		}
		return term1;
	}

	int apb = a + b;
	int term2 = logexp_device((apb > 0 ? apb : -apb), Dint1, Dint2, Dint3);
	int amb = a - b;
	int term3 = logexp_device((amb > 0 ? amb : -amb), Dint1, Dint2, Dint3);
	int result = term1 + term2 - term3;

	// Don't abort when overflowing, just saturate the QLLR
	if (result > QLLR_MAX) {
		return QLLR_MAX;
	}
	if (result < -QLLR_MAX) {
		return -QLLR_MAX;
	}
	return result;
}

__global__ 
void updateCheckNode_kernel( const int ncheck, const int nvar, 
	const int* sumX2, const int* mvc, const int* jind, 
	const short int Dint1, const short int Dint2, const short int Dint3, 
	int* d_ml, int* d_mr, const int max_cnd, const int QLLR_MAX,
	int* mcv )
{	//	mvc const(input)-> mcv (output)
	int j = blockIdx.x * blockDim.x + threadIdx.x;

	if( j>= ncheck )
		return;

	int ml[MAX_CHECK_NODE];//int* ml	= d_ml	+ j * max_cnd;
	int mr[MAX_CHECK_NODE];//int* mr	= d_mr	+ j * max_cnd;

	int nodes = sumX2[j];

	nodes--;

	// compute partial sums from the left and from the right
#if USE_TEXTURE_ADDRESS
	int indexL = jind[j];
	ml[0] = tex2D(texMVC, indexL%nvar, indexL/nvar);
	int indexR = jind[j+nodes*ncheck];
	mr[0] = tex2D(texMVC, indexR%nvar, indexR/nvar);
	for(int i = 1; i < nodes; i++ ) {
		indexL = jind[j+i*ncheck];
		int valueL = tex2D(texMVC, indexL%nvar, indexL/nvar);;
		indexR = jind[j+(nodes-i)*ncheck];
		int valueR = tex2D(texMVC, indexR%nvar, indexR/nvar);;

		ml[i] = Boxplus( ml[i-1], valueL, Dint1, Dint2, Dint3, QLLR_MAX );
		mr[i] = Boxplus( mr[i-1], valueR, Dint1, Dint2, Dint3, QLLR_MAX );
	}
#else
	ml[0] = mvc[jind[j]];
	mr[0] = mvc[jind[j+nodes*ncheck]];
	for(int i = 1; i < nodes; i++ ) {
		ml[i] = Boxplus( ml[i-1], mvc[jind[j+i*ncheck]], Dint1, Dint2, Dint3, QLLR_MAX );
		mr[i] = Boxplus( mr[i-1], mvc[jind[j+(nodes-i)*ncheck]], Dint1, Dint2, Dint3, QLLR_MAX );
	}
#endif
	// merge partial sums
	mcv[j] = mr[nodes-1];
	mcv[j+nodes*ncheck] = ml[nodes-1];
	for(int i = 1; i < nodes; i++ )
		mcv[j+i*ncheck] = Boxplus( ml[i-1], mr[nodes-1-i], Dint1, Dint2, Dint3, QLLR_MAX );

}

__global__ 
void updateCheckNodeShared_kernel( const int ncheck, const int nvar, 
	const int* sumX2, const int* mvc, const int* jind, 
	const short int Dint1, const short int Dint2, const short int Dint3, 
	int* d_ml, int* d_mr, const int max_cnd, const int QLLR_MAX,
	int* mcv )
{	//	mvc const(input)-> mcv (output)
	int j = blockIdx.x * blockDim.x + threadIdx.x;

	if( j>= ncheck )
		return;
#if USE_SHARED_MLR
	__shared__ int s_ml[MAX_CHECK_NODE*SIZE_BLOCK];
	int* ml	= s_ml	+ threadIdx.x * MAX_CHECK_NODE;
	__shared__ int s_mr[MAX_CHECK_NODE*SIZE_BLOCK];
	int* mr	= s_mr	+ threadIdx.x * MAX_CHECK_NODE;
#else
	int ml[MAX_CHECK_NODE];
	int mr[MAX_CHECK_NODE];
#endif

	int nodes = sumX2[j];

	nodes--;

	// compute partial sums from the left and from the right
#if USE_SHARED_MVC
	__shared__ int s_mvc[MAX_CHECK_NODE*SIZE_BLOCK];
	for(int i = 0; i < sumX2[j]; i++ ) {
		s_mvc[threadIdx.x+i*SIZE_BLOCK] = mvc[ jind[j+i*ncheck] ];
	}

	ml[0] = s_mvc[threadIdx.x];
	mr[0] = s_mvc[threadIdx.x+nodes*SIZE_BLOCK];
	for(int i = 1; i < nodes; i++ ) {
		ml[i] = Boxplus( ml[i-1], s_mvc[threadIdx.x+i*SIZE_BLOCK], Dint1, Dint2, Dint3, QLLR_MAX );
		mr[i] = Boxplus( mr[i-1], s_mvc[threadIdx.x+(nodes-i)*SIZE_BLOCK], Dint1, Dint2, Dint3, QLLR_MAX );
	}
#else
	ml[0] = mvc[jind[j]];
	mr[0] = mvc[jind[j+nodes*ncheck]];
	for(int i = 1; i < nodes; i++ ) {
		ml[i] = Boxplus( ml[i-1], mvc[jind[j+i*ncheck]], Dint1, Dint2, Dint3, QLLR_MAX );
		mr[i] = Boxplus( mr[i-1], mvc[jind[j+(nodes-i)*ncheck]], Dint1, Dint2, Dint3, QLLR_MAX );
	}
#endif

	// merge partial sums
	mcv[j] = mr[nodes-1];
	mcv[j+nodes*ncheck] = ml[nodes-1];
	for(int i = 1; i < nodes; i++ )
		mcv[j+i*ncheck] = Boxplus( ml[i-1], mr[nodes-1-i], Dint1, Dint2, Dint3, QLLR_MAX );

}

__global__ 
void initializeMVC_kernel(const int nvar, 
	const int* sumX1,
	const int* LLRin,
	int* mvc) 
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	if( i>= nvar )
		return;

	int index = i;
    for (int j = 0; j < sumX1[i]; j++) {
      mvc[index] = LLRin[i];
      index += nvar;
    }
}


__global__ 
void updateVariableNodeAndCheck_kernel( const int nvar, const int ncheck, 
	const int* sumX1, const int* sumX2, const int* iind, const int* d_V, 
	const int * LLRin, const int* mcv, 
	int * LLRout, int* mvc,
	int* d_synd ) 
{	//	mcv const(input)-> mvc (output)
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	
	if( i>= nvar )
		return;

	int mvc_temp = LLRin[i];
	for (int jp = 0; jp < sumX1[i]; jp++) {
		int index = iind[i + jp*nvar];
		mvc_temp +=  mcv[index];
	}
	LLRout[i] = mvc_temp;
	int index_iind = i;  // tracks i+j*nvar
	for (int jj = 0; jj < sumX1[i]; jj++) {
		mvc[index_iind] = mvc_temp - mcv[iind[index_iind]];
		index_iind += nvar;
	}
	

	__syncthreads();

	int j = i;

	if( j>= ncheck )
		return;

	if( j == 0 )
		*d_synd = 0;

	__syncthreads();

	int ii, vi;
	int synd = 0;
	int vind = j; // tracks j+i*ncheck
	for (ii = 0; ii < sumX2[j]; ii++) {
		vi = d_V[vind];
		if (LLRout[vi] < 0) {
			synd++;
		}
		vind += ncheck;
	}

	if( synd&1 )
		atomicAdd( d_synd, 1 );	// synd is even ?
}