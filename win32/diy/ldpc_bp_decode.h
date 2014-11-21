
#pragma once


bool syndrome_check(int* LLR,
	int ncheck, 
	int* sumX2, 
	int* V) ;

int  logexp(int x,
	short int Dint1, short int Dint2, short int Dint3,	//! Decoder (lookup-table) parameters
	int* logexp_table );		//! The lookup tables for the decoder

int Boxplus(int a, int b,
	short int Dint1, short int Dint2, short int Dint3,	//! Decoder (lookup-table) parameters
	int* logexp_table );		//! The lookup tables for the decoder

int bp_decode(int *LLRin, int *LLRout,
	int nvar, int ncheck, 
	int* V, int* sumX1, int* sumX2, int* iind, int* jind,	// Parity check matrix parameterization
	int* mvc, int* mcv,	// temporary storage for decoder (memory allocated when codec defined)
	//LLR_calc_unit& llrcalc,		//!< LLR calculation unit
	short int Dint1, short int Dint2, short int Dint3,	//! Decoder (lookup-table) parameters
	int* logexp_table,		//! The lookup tables for the decoder
	bool psc = true,			//!< check syndrom after each iteration
	int max_iters = 50 );		//!< Maximum number of iterations
