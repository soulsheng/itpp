
#pragma	once


class	driverUpdataVar
{

public:
	bool	launch();		//!	launch kernel 

	bool	verify();	//! verify whether the result is right or wrong
	
	driverUpdataVar();

	~driverUpdataVar();

private:
	int nvar, ncheck;
	int nmaxX1, nmaxX2;

	// host
	int *sumX1, *sumX2;
	int *iind, *jind;
	int *mvc, *mcv;
	int *input;	char *output;

	int *ref_mvc;	char *ref_output;

	// device
	int *d_sumX1, *d_sumX2;
	int *d_iind, *d_jind;
	int *d_mvc, *d_mcv;
	int *d_input;	char *d_output;

};
