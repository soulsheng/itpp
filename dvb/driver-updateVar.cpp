
#include "driver-updateVar.cuh"

#include <iostream>
using	namespace	std;


void main()
{
	bool bStatus = false;

	driverUpdataVar ldpc;

	if( !ldpc.launch() )
		cout << "failed to launch" << endl;

	if( !ldpc.verify() )
		cout << "failed to verify" << endl;

}

