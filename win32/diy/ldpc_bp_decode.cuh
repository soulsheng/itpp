
#pragma once

class ldpc_gpu
{
protected:
bool syndrome_check_gpu( ) ;

void updateVariableNode_gpu( ) ;

void updateCheckNode_gpu( );

void initializeMVC_gpu( );

public:
int bp_decode(int *LLRin, int *LLRout,
	bool psc = true,			//!< check syndrom after each iteration
	int max_iters = 50 );		//!< Maximum number of iterations

int bp_decode_once(int *LLRin, int *LLRout,
	bool psc = true,			//!< check syndrom after each iteration
	int max_iters = 50 );		//!< Maximum number of iterations

	bool	initialize( int nvar, int ncheck,
	int nmaxX1, int nmaxX2,
	int* V, int* sumX1, int* sumX2, int* iind, int* jind,	// Parity check matrix parameterization
	int* mvc, int* mcv,	// temporary storage for decoder (memory allocated when codec defined)
	short int Dint1, short int Dint2, short int Dint3,
	int* logexp_table		//! The lookup tables for the decoder
	);

	~ldpc_gpu();

private:
	bool	release();

private:
	int* d_synd ;

	int* d_sumX1 ;
	int* d_sumX2 ;
	int* d_mcv ;
	int* d_mvc ;
	int* d_iind ;
	int* d_jind ;
	int* d_V ;

	//int* d_logexp_table ;
	
	int *d_ml, *d_mr ;
	
	int* d_LLRin ;
	int* d_LLRout ;
	
	char*	d_bLLR;	// 1 (LLR < 0),  0 (LLR >= 0)

private:
	int nvar, ncheck;
	int nmaxX1, nmaxX2; // max(sumX1) max(sumX2)
	short int Dint1, Dint2, Dint3;	//! Decoder (lookup-table) parameters
	int max_cnd;	//! Maximum check node degree that the class can handle
	int QLLR_MAX;
};