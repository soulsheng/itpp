
#pragma	once


class	driverUpdataChk
{

public:
	bool	launch();		//!	launch kernel 

	bool	verify();	//! verify whether the result is right or wrong
	
	driverUpdataChk();

	~driverUpdataChk();

private:
	int nvar, ncheck;
	int nmaxX1, nmaxX2;
	short int Dint1, Dint2, Dint3;	//! Decoder (lookup-table) parameters
	int QLLR_MAX;

	// host
	int *sumX1, *sumX2;
	int *iind, *jind;
	int *mvc, *mcv;
	int *input;	char *output;

	int *ref_mcv;	char *ref_output;

	int *logexp_table;

	// device
	int *d_sumX1, *d_sumX2;
	int *d_iind, *d_jind;
	int *d_mvc, *d_mcv;
	int *d_input;	char *d_output;

	int* d_logexp_table ;

};
