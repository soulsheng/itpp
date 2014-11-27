
#pragma once

#define		TABLE_SIZE_DINT2	300
#define		MAX_CHECK_NODE		100
#define		TABLE_SIZE_CODE		16200
#define		USE_TABLE_CODE		0
#define		USE_TEXTURE_ADDRESS	0

#if USE_TEXTURE_ADDRESS
texture<int, 2, cudaReadModeElementType> texMCV;
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
	
	atomicAdd( d_synd,  synd&1 );	// synd is even ?
}

__global__ 
void updateVariableNode_kernel( const int nvar, const int ncheck, const int* sumX1, const int* mcv, const int* iind, const int * LLRin, 
	int * LLRout, int* mvc, char* bLLR ) 
{	//	mcv const(input)-> mvc (output)
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	
	if( i>= nvar )
		return;
#if 0
	switch (sumX1[i]) {
		case 0:
			return;// cout << "LDPC_Code::bp_decode(): sumX1[i]=0" << endl;
		case 1: {
			/* This case is rare but apparently occurs for codes used in
			the DVB-T2 standard.
			*/
			int m0 = mcv[iind[i]];
			mvc[i] = LLRin[i];
			LLRout[i] = LLRin[i] + m0;
			break;
				}
		case 2: {
			int m0 = mcv[iind[i]];
			int i1 = i + nvar;
			int m1 = mcv[iind[i1]];
			mvc[i] = LLRin[i] + m1;
			mvc[i1] = LLRin[i] + m0;
			LLRout[i] = mvc[i1] + m1;
			break;
				}
		case 3: {
			int i0 = i;
			int m0 = mcv[iind[i0]];
			int i1 = i0 + nvar;
			int m1 = mcv[iind[i1]];
			int i2 = i1 + nvar;
			int m2 = mcv[iind[i2]];
			LLRout[i] = LLRin[i] + m0 + m1 + m2;
			mvc[i0] = LLRout[i] - m0;
			mvc[i1] = LLRout[i] - m1;
			mvc[i2] = LLRout[i] - m2;
			break;
				}
		case 4: {
			int i0 = i;
			int m0 = mcv[iind[i0]];
			int i1 = i0 + nvar;
			int m1 = mcv[iind[i1]];
			int i2 = i1 + nvar;
			int m2 = mcv[iind[i2]];
			int i3 = i2 + nvar;
			int m3 = mcv[iind[i3]];
			LLRout[i] = LLRin[i] + m0 + m1 + m2 + m3;
			mvc[i0] = LLRout[i] - m0;
			mvc[i1] = LLRout[i] - m1;
			mvc[i2] = LLRout[i] - m2;
			mvc[i3] = LLRout[i] - m3;
			break;
				}
		default:   { // differential update
#endif
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
		//		   }
		//}
	
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
void updateCheckNode_kernel( const int ncheck, 
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
#if 0
	switch (sumX2[j]) {
		case 0:
			return;//cout << "LDPC_Code::bp_decode(): sumX2[j]=0" << endl;
		case 1:
			return;//cout << "LDPC_Code::bp_decode(): sumX2[j]=1" << endl;
		case 2: {
			mcv[j+ncheck] = mvc[jind[j]];
			mcv[j] = mvc[jind[j+ncheck]];
			break;
				}
		case 3: {
			int j0 = j;
			int m0 = mvc[jind[j0]];
			int j1 = j0 + ncheck;
			int m1 = mvc[jind[j1]];
			int j2 = j1 + ncheck;
			int m2 = mvc[jind[j2]];
			mcv[j0] = Boxplus(m1, m2, Dint1, Dint2, Dint3, QLLR_MAX);
			mcv[j1] = Boxplus(m0, m2, Dint1, Dint2, Dint3, QLLR_MAX);
			mcv[j2] = Boxplus(m0, m1, Dint1, Dint2, Dint3, QLLR_MAX);
			break;
				}
		case 4: {
			int j0 = j;
			int m0 = mvc[jind[j0]];
			int j1 = j0 + ncheck;
			int m1 = mvc[jind[j1]];
			int j2 = j1 + ncheck;
			int m2 = mvc[jind[j2]];
			int j3 = j2 + ncheck;
			int m3 = mvc[jind[j3]];
			int m01 = Boxplus(m0, m1, Dint1, Dint2, Dint3, QLLR_MAX);
			int m23 = Boxplus(m2, m3, Dint1, Dint2, Dint3, QLLR_MAX);
			mcv[j0] = Boxplus(m1, m23, Dint1, Dint2, Dint3, QLLR_MAX);
			mcv[j1] = Boxplus(m0, m23, Dint1, Dint2, Dint3, QLLR_MAX);
			mcv[j2] = Boxplus(m01, m3, Dint1, Dint2, Dint3, QLLR_MAX);
			mcv[j3] = Boxplus(m01, m2, Dint1, Dint2, Dint3, QLLR_MAX);
			break;
				}
		case 5: {
			int j0 = j;
			int m0 = mvc[jind[j0]];
			int j1 = j0 + ncheck;
			int m1 = mvc[jind[j1]];
			int j2 = j1 + ncheck;
			int m2 = mvc[jind[j2]];
			int j3 = j2 + ncheck;
			int m3 = mvc[jind[j3]];
			int j4 = j3 + ncheck;
			int m4 = mvc[jind[j4]];
			int m01 = Boxplus(m0, m1, Dint1, Dint2, Dint3, QLLR_MAX);
			int m02 = Boxplus(m01, m2, Dint1, Dint2, Dint3, QLLR_MAX);
			int m34 = Boxplus(m3, m4, Dint1, Dint2, Dint3, QLLR_MAX);
			int m24 = Boxplus(m2, m34, Dint1, Dint2, Dint3, QLLR_MAX);
			mcv[j0] = Boxplus(m1, m24, Dint1, Dint2, Dint3, QLLR_MAX);
			mcv[j1] = Boxplus(m0, m24, Dint1, Dint2, Dint3, QLLR_MAX);
			mcv[j2] = Boxplus(m01, m34, Dint1, Dint2, Dint3, QLLR_MAX);
			mcv[j3] = Boxplus(m02, m4, Dint1, Dint2, Dint3, QLLR_MAX);
			mcv[j4] = Boxplus(m02, m3, Dint1, Dint2, Dint3, QLLR_MAX);
			break;
				}
		case 6: {
			int j0 = j;
			int m0 = mvc[jind[j0]];
			int j1 = j0 + ncheck;
			int m1 = mvc[jind[j1]];
			int j2 = j1 + ncheck;
			int m2 = mvc[jind[j2]];
			int j3 = j2 + ncheck;
			int m3 = mvc[jind[j3]];
			int j4 = j3 + ncheck;
			int m4 = mvc[jind[j4]];
			int j5 = j4 + ncheck;
			int m5 = mvc[jind[j5]];
			int m01 = Boxplus(m0, m1, Dint1, Dint2, Dint3, QLLR_MAX);
			int m23 = Boxplus(m2, m3, Dint1, Dint2, Dint3, QLLR_MAX);
			int m45 = Boxplus(m4, m5, Dint1, Dint2, Dint3, QLLR_MAX);
			int m03 = Boxplus(m01, m23, Dint1, Dint2, Dint3, QLLR_MAX);
			int m25 = Boxplus(m23, m45, Dint1, Dint2, Dint3, QLLR_MAX);
			int m0145 = Boxplus(m01, m45, Dint1, Dint2, Dint3, QLLR_MAX);
			mcv[j0] = Boxplus(m1, m25, Dint1, Dint2, Dint3, QLLR_MAX);
			mcv[j1] = Boxplus(m0, m25, Dint1, Dint2, Dint3, QLLR_MAX);
			mcv[j2] = Boxplus(m0145, m3, Dint1, Dint2, Dint3, QLLR_MAX);
			mcv[j3] = Boxplus(m0145, m2, Dint1, Dint2, Dint3, QLLR_MAX);
			mcv[j4] = Boxplus(m03, m5, Dint1, Dint2, Dint3, QLLR_MAX);
			mcv[j5] = Boxplus(m03, m4, Dint1, Dint2, Dint3, QLLR_MAX);
			break;
				}
		default: {
#endif
			int nodes = sumX2[j];

			nodes--;

			// compute partial sums from the left and from the right
			ml[0] = mvc[jind[j]];
			mr[0] = mvc[jind[j+nodes*ncheck]];
			for(int i = 1; i < nodes; i++ ) {
				ml[i] = Boxplus( ml[i-1], mvc[jind[j+i*ncheck]], Dint1, Dint2, Dint3, QLLR_MAX );
				mr[i] = Boxplus( mr[i-1], mvc[jind[j+(nodes-i)*ncheck]], Dint1, Dint2, Dint3, QLLR_MAX );
			}

			// merge partial sums
			mcv[j] = mr[nodes-1];
			mcv[j+nodes*ncheck] = ml[nodes-1];
			for(int i = 1; i < nodes; i++ )
				mcv[j+i*ncheck] = Boxplus( ml[i-1], mr[nodes-1-i], Dint1, Dint2, Dint3, QLLR_MAX );
				 //}
		//}  // switch statement

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