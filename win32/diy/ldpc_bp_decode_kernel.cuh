
#pragma once

__global__ 
void syndrome_check_kernel(int *d_LLR,
	int* d_sumX2, int ncheck, 
	int* d_V,
	int* d_synd) 
{
	int j = blockIdx.x * blockDim.x + threadIdx.x;

	if( j>= ncheck )
		return;

	int i, vi;

	int vind = j; // tracks j+i*ncheck
	for (i = 0; i < d_sumX2[j]; i++) {
		vi = d_V[vind];
		if (d_LLR[vi] < 0) {
			d_synd[j]++;
		}
		vind += ncheck;
	}
	
	d_synd[j] = !(d_synd[j]&1);	// d_synd[j] is even ?
}

__global__ 
void updateVariableNode_kernel( int nvar, int* sumX1, int* mcv, int* mvc, int* iind, int * LLRin, int * LLRout ) 
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	
	if( i>= nvar )
		return;

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
			int mvc_temp = LLRin[i];
			int index_iind = i; // tracks i+jp*nvar
			for (int jp = 0; jp < sumX1[i]; jp++) {
				mvc_temp +=  mcv[iind[index_iind]];
				index_iind += nvar;
			}
			LLRout[i] = mvc_temp;
			index_iind = i;  // tracks i+j*nvar
			for (int j = 0; j < sumX1[i]; j++) {
				mvc[index_iind] = mvc_temp - mcv[iind[index_iind]];
				index_iind += nvar;
			}
				   }
		}
}