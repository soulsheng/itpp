
#pragma once

//! Maximum value of vector
int max(int *v, int N);
//! Minimum value of vector
int min(int *v, int N);

class ldpc_decoder{
public:

bool syndrome_check(char* LLR,
	int ncheck, 
	int* sumX2, 
	int* V) ;

int  logexp(int x,
	short int Dint1, short int Dint2, short int Dint3,	//! Decoder (lookup-table) parameters
	int* logexp_table );		//! The lookup tables for the decoder

int Boxplus(int a, int b,
	short int Dint1, short int Dint2, short int Dint3,	//! Decoder (lookup-table) parameters
	int* logexp_table );		//! The lookup tables for the decoder

void initialize(	int nvar, int ncheck, 
	int nmaxX1, int nmaxX2, // max(sumX1) max(sumX2)
	int* V, int* sumX1, int* sumX2, int* iind, int* jind,	// Parity check matrix parameterization
	int* mvc, int* mcv,	// temporary storage for decoder (memory allocated when codec defined)
	//LLR_calc_unit& llrcalc,		//!< LLR calculation unit
	short int Dint1, short int Dint2, short int Dint3,	//! Decoder (lookup-table) parameters
	int* logexp_table,		//! The lookup tables for the decoder
	bool psc = true,			//!< check syndrom after each iteration
	int max_iters = 50 );		//!< Maximum number of iterations

void updateCheckNode( int ncheck, int* sumX2, 
	int* mcv, int* mvc, int* jind, 
	short int Dint1, short int Dint2, short int Dint3, 
	int* logexp_table );

void updateVariableNode( int nvar, int* sumX1, 
	int* mcv, int* mvc, int* iind, 
	int * LLRin, char * LLRout );

void initializeMVC( int nvar, int* sumX1, int* mvc, int * LLRin );

int bp_decode(int *LLRin, char *LLRout);

protected:
	int *LLRin; char *LLRout;
	int nvar, ncheck;
	int nmaxX1, nmaxX2; // max(sumX1) max(sumX2)
	int* V, * sumX1, * sumX2, * iind, * jind;	// Parity check matrix parameterization
	int* mvc; int* mcv;	// temporary storage for decoder (memory allocated when codec defined)
	short int Dint1, Dint2, Dint3;	//! Decoder (lookup-table) parameters
	int* logexp_table;		//! The lookup tables for the decoder
	bool psc;			//!< check syndrom after each iteration
	int max_iters;
};