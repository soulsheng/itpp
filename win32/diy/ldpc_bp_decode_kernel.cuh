
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

