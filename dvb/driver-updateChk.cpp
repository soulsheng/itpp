
#include "driver-updateChk.cuh"

#include <cuda_runtime.h>
#include <iostream>
using	namespace	std;


void main()
{
	bool bStatus = false;

	driverUpdataChk ldpc;

	if( !ldpc.launch() )
		cout << "Failed to launch" << endl;

	if( !ldpc.verify() )
		cout << "Failed to verify" << endl;
	else
		cout << "Succeed to launch cuda kernel and verify that result is right" << endl;

	cudaDeviceReset();
	system( "pause" );

}

