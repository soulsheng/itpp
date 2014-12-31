
#include "driver-updateChk.cuh"

#include <cuda_runtime.h>
#include <iostream>
using	namespace	std;


void main()
{
	bool bStatus = false;

	driverUpdataChk ldpc;

	if( !ldpc.launch() )
		cout << "failed to launch" << endl;

	if( !ldpc.verify() )
		cout << "failed to verify" << endl;
	else
		cout << "succeed to launch and verify that result is right" << endl;

	cudaDeviceReset();
	system( "pause" );

}

