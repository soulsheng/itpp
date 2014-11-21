
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


int bp_decode_gpu(int *LLRin, int *LLRout,
	int nvar, int ncheck, 
	int nmaxX1, int nmaxX2, // max(sumX1) max(sumX2)
	int* V, int* sumX1, int* sumX2, int* iind, int* jind,	// Parity check matrix parameterization
	int* mvc, int* mcv,	// temporary storage for decoder (memory allocated when codec defined)
	//LLR_calc_unit& llrcalc,		//!< LLR calculation unit
	short int Dint1, short int Dint2, short int Dint3,	//! Decoder (lookup-table) parameters
	int* logexp_table,		//! The lookup tables for the decoder
	bool psc = true,			//!< check syndrom after each iteration
	int max_iters = 50 );		//!< Maximum number of iterations
