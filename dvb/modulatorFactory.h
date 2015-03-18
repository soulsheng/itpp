
#pragma once

#include <itpp/itcomm.h>
using namespace itpp;

#include <vector>
#include <map>
#include <iostream>
using namespace std;

enum	MOD_TYPE
{
	MOD_DEFAULT,//	0
	MOD_BPSK,	//	1
	MOD_QPSK,	//	2
	MOD_8PSK,	//	3
	MOD_16APSK,	//	4
	MOD_32APSK,	//	5
	MOD_COUNT	//	6
};

enum CODE_RATE
{
	C2_3,
	C3_4,
	C4_5,
	C5_6,
	C8_9,
	C9_10
};

enum FRAME_TYPE
{
	FECFRAME_NORMAL,
	FECFRAME_SHORT
};

struct SymbolTable
{
	cvec	symbols;
	ivec	bits10symbols;
	Vec<bvec>	bits2symbols;
	int		k;
	double  r1, r2, r3;
	CODE_RATE rate;
	FRAME_TYPE framesize;

	SymbolTable(int k, CODE_RATE rate=C3_4, FRAME_TYPE framesize=FECFRAME_NORMAL);
	cvec&	getSymbols()	{ return symbols;}
	ivec&	getBits10Symbols()	{ return bits10symbols;}
};


class ModulatorFactory
{
public:
	ModulatorFactory();
	~ModulatorFactory();

	Modulator_2D* findModulator(MOD_TYPE modType);

protected:
private:
	typedef map<MOD_TYPE, Modulator_2D*> ModPool;
	typedef map<MOD_TYPE, Modulator_2D*>::iterator ModPoolItr;
	typedef pair<MOD_TYPE, Modulator_2D*> ModPoolPair;

	ModPool	m_modPool;
};
