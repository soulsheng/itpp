
#pragma once


bool syndrome_check_gpu(int* LLR, int nvar,  
	int* sumX2, int ncheck,
	int* V, int nmaxX2) ;

void updateVariableNode_gpu( int nvar, int ncheck, int nmaxX1, int nmaxX2, 
	int* sumX1, int* mcv, int* mvc, int* iind, int * LLRin, int * LLRout ) ;
