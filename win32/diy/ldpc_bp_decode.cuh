
#pragma once


bool syndrome_check_gpu(int* LLR, int nvar,  
	int* sumX2, int ncheck,
	int* V, int nmaxX2) ;

void updateVariableNode_gpu( int nvar, int ncheck, int nmaxX1, int nmaxX2, 
	int* sumX1, int* mcv, int* mvc, int* iind, int * LLRin, int * LLRout ) ;

void updateCheckNode_gpu( int nvar, int ncheck, int nmaxX1, int nmaxX2, 
	int* sumX2, int* mcv, int* mvc, int* jind, 
	short int Dint1, short int Dint2, short int Dint3, int* logexp_table,
	int* jj, int* m, int* ml, int* mr, int max_cnd, int QLLR_MAX );