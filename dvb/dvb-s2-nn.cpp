#include <itpp/itcomm.h>
#include <sstream>
#include <cuda_runtime.h>
#include "ldpc_bp_decode.h"
#include "ldpc_bp_decode.cuh"
#include "helper_timer.h"
//#include "driverUtility.h"
#include "dvbUtility.h"
#include "modulatorFactory.h"
#include "bbheader_bb.h"
#include "bbscrambler_bb.h"
#include "bch_bm.h"

using namespace std;
using namespace itpp;


#define		TIME_STEP		3	

#define		USE_GPU			0
#define		USE_ALIST		0


int main(int argc, char **argv)
{
#if 1
	FILE	*fp;
	fp = fopen( FILENAME_CODE_DOWNLOAD, "rb" );
	if( fp == NULL )
	{
		cout << FILENAME_CODE_DOWNLOAD << " can not open " << endl;
		return 0;
	}
#else
	ifstream  bitDownload;
	bitDownload.open( FILENAME_CODE_DOWNLOAD );
#endif

	CODE_RATE	rate = C3_4;
	FRAME_TYPE	type = FECFRAME_NORMAL;
	// step 0: intialize ldpc,bch,bpsk,awgn
#if USE_ALIST
	LDPC_Parity	lp( FILENAME_ALIST, "alist" );
	LDPC_Code	ldpc( &lp );
#else
	LDPC_Generator_Systematic G; // for codes created with ldpc_gen_codes since generator exists
	LDPC_Code ldpc(FILENAME_IT34, &G);
#endif

	BCH_BM	bch;
	bch.initialize();
	bch.setCode( C3_4, FECFRAME_NORMAL );
	double N0 = pow(10.0, -EBNO / 10.0) / ldpc.get_rate();

	ModulatorFactory	mods;	// 调制解调器件库

	BERC berc;  // Counters for coded and uncoded BER
	BLERC per(SIZE_PACKET); // Counter for coded FER

	RNG_randomize();

	int nldpc = ldpc.get_nvar();
	int kldpc = ldpc.get_ninfo();             // number of bits per codeword

	int nSplit = kldpc / N_BCH;

	//int Kbch = nSplit * K_BCH;

	// print parameter value
	int nmaxX1 = max(ldpc.sumX1._data(), ldpc.sumX1.size());
	int nmaxX2 = max(ldpc.sumX2._data(), ldpc.sumX2.size());
	int nminX1 = min(ldpc.sumX1._data(), ldpc.sumX1.size());
	int nminX2 = min(ldpc.sumX2._data(), ldpc.sumX2.size());

	int nmaxI = max(ldpc.iind._data(), ldpc.iind.size());
	int nmaxJ = max(ldpc.jind._data(), ldpc.jind.size());
	int nminI = min(ldpc.iind._data(), ldpc.iind.size());
	int nminJ = min(ldpc.jind._data(), ldpc.jind.size());

#if 1
	cout << "max(iind) = " << nmaxI << endl;// max(iind) = nvar*nmaxX1-1
	cout << "max(jind) = " << nmaxJ << endl;// max(jind) = nvar*nmaxX1-1
	cout << "min(iind) = " << nminI << endl;// min(iind) = 0
	cout << "min(jind) = " << nminJ << endl;// min(jind) = 0

	cout << "ldpc.nvar = " << ldpc.nvar << endl;		// nvar = 16200
	cout << "ldpc.ncheck = " << ldpc.ncheck << endl;	// ncheck = 8100//8073 
	cout << "ldpc.sumX1.size() = " << ldpc.sumX1.size() << endl;	// = nvar
	cout << "ldpc.sumX2.size() = " << ldpc.sumX2.size() << endl;	// = ncheck
	cout << "max(sumX1) = " << nmaxX1 << endl;// max(sumX1) = 3//19
	cout << "max(sumX2) = " << nmaxX2 << endl;// max(sumX2) = 6//10
	cout << "min(sumX1) = " << nminX1 << endl;// min(sumX1) = 3//2
	cout << "min(sumX2) = " << nminX2 << endl;// min(sumX2) = 6//7
	cout << "ldpc.V.size() = " << ldpc.V.size() << endl;			// = ncheck * max(sumX2)
	cout << "ldpc.iind.size() = " << ldpc.iind.size() << endl;		// = nvar * max(sumX1)
	cout << "ldpc.jind.size() = " << ldpc.jind.size() << endl;		// = ncheck * max(sumX2)

	cout << "ldpc.mvc.size() = " << ldpc.mvc.size() << endl;		// = nvar * max(sumX1)
	cout << "ldpc.mcv.size() = " << ldpc.mcv.size() << endl;		// = ncheck * max(sumX2)

	cout << "ldpc.llrcalc.Dint1 = " << ldpc.llrcalc.Dint1 << endl;	// Dint1 = 12
	cout << "ldpc.llrcalc.Dint2 = " << ldpc.llrcalc.Dint2 << endl;	// Dint2 = 300
	cout << "ldpc.llrcalc.Dint3 = " << ldpc.llrcalc.Dint3 << endl;	// Dint3 = 7

	cout << "ldpc.llrcalc.logexp_table.size() = " << ldpc.llrcalc.logexp_table.size() << endl;// = 300
#endif

#if		USE_GPU
	ldpc_gpu	ldpc_gpu_diy;
	ldpc_gpu_diy.initialize(ldpc.nvar, ldpc.ncheck, 
		nmaxX1, nmaxX2, 
		ldpc.sumX1._data(), ldpc.sumX2._data(), ldpc.iind._data(), ldpc.jind._data(), ldpc.V._data(), 	// Parity check matrix parameterization
		ldpc.llrcalc.Dint1, ldpc.llrcalc.Dint2, ldpc.llrcalc.Dint3,	//! Decoder (lookup-table) parameters
		ldpc.llrcalc.logexp_table._data());
#endif

	char bitOut[FRAME_SIZE_NORMAL]={0};
	char bitOutTmp[FRAME_SIZE_NORMAL]={0};


	int COUNT_REPEAT = COUNT_REPEAT_DEF;

	char *bitsPacketsPadding = new char[nldpc];
	char *bitsLDPC = new char[nldpc];
	double *softbits_buf = new double[nldpc];
	char messageRef[FRAME_SIZE_NORMAL]={0};		// information received

	vec			timerValue(COUNT_REPEAT);

	vec			timerStepValue(COUNT_REPEAT*TIME_STEP);
	ivec		countIteration(COUNT_REPEAT);



	cout << "Demodulating and decoding the bit stream !" << endl ;
	cout << "Please wait for a few seconds ..." << endl << endl;

	// header
	bbheader_bb*	pBBHeader = bbheader_bb::make(C1_2, RO_0_20, FECFRAME_SHORT);
	bbscrambler_bb*	pBBScrambler = bbscrambler_bb::make(C1_2, FECFRAME_SHORT);

	int Kbch = bch.getK();

	for (int64_t i = 0; i < COUNT_REPEAT; i ++) 
	{
		// demodulate

		int nModTypeRAND = MOD_16APSK;
#if 0
		bitfileMOD >> nModTypeRAND ;
		cout << "nModTypeRAND = " << nModTypeRAND << endl;
#endif
		MOD_TYPE	modType = (MOD_TYPE)nModTypeRAND;
		Modulator_2D* pModulator = mods.findModulator( modType );

		if ( NULL == pModulator )
		{
			cout << "can not modulate" << endl;
			continue;
		}

		int	nSizeMod = nldpc*2/pModulator->get_k();
		double *bitsMOD  = new double[nSizeMod];
		short	*bitsMODS  = new short[nSizeMod];

		fread( bitsMODS, sizeof(short), nSizeMod, fp );
		for ( int j=0; j<nSizeMod; j++ )
		{
			bitsMOD[j] = bitsMODS[j]/float(1<<13) ;
		}

		cvec	cAWGN( nSizeMod/2 );
		convertBufferToVec( bitsMOD, cAWGN );

		cout << "cAWGN.left( 10 )" << cAWGN.left( 10 ) << endl;

		vec softbits = pModulator->demodulate_soft_bits(cAWGN, N0);

		// ldpc Decode the received bits
		QLLRvec llrIn = ldpc.get_llrcalc().to_qllr(softbits);


#if		USE_GPU
		countIteration[i] = ldpc_gpu_diy.bp_decode_once( llrIn._data(), bitOut ); 
#else
		countIteration[i] = bp_decode( llrIn._data(), bitOut, 
			ldpc.nvar, ldpc.ncheck, 
			nmaxX1, nmaxX2, 
			ldpc.V._data(), ldpc.sumX1._data(), ldpc.sumX2._data(), ldpc.iind._data(), ldpc.jind._data(),	// Parity check matrix parameterization
			ldpc.mvc._data(), ldpc.mcv._data(),	// temporary storage for decoder (memory allocated when codec defined)
			//ldpc.llrcalc );		//!< LLR calculation unit
			ldpc.llrcalc.Dint1, ldpc.llrcalc.Dint2, ldpc.llrcalc.Dint3,	//! Decoder (lookup-table) parameters
			ldpc.llrcalc.logexp_table._data());		//! The lookup tables for the decoder

#endif

		bvec bitsoutLDPCDec(nldpc);
		for (int j=0;j<nldpc;j++)
			bitsoutLDPCDec[j] = bitOut[j];

		pBBScrambler->work( Kbch, bitOut, bitOutTmp );

		bool bStatus = pBBHeader->verify( bitOutTmp );

		bvec bitsoutDec(Kbch-80);
		convertBufferToVec( bitOutTmp+80 , bitsoutDec );

		// step 9: verify result, Count the number of errors
		bvec bitsinEnc(Kbch-80);
		convertBufferToVec( messageRef, bitsinEnc );

		cout << "bitsoutLDPCDec.left(100)" << bitsoutLDPCDec.left(100) << endl << endl;
		cout << "bitsoutDec.left(100)" << bitsoutDec.left(100) << endl << endl;

		berc.count(bitsinEnc, bitsoutDec);
		per.count(bitsinEnc, bitsoutDec);

		free( bitsMOD );

    }

	cout << "Eb/N0 = " << EBNO << endl << "  Simulated "
		<< COUNT_REPEAT << " frames made of "
		<< berc.get_total_bits() << " bits. "
		<< per.get_total_blocks() << " packets. " << endl
		<< "Obtained " << berc.get_errors() << " bit errors. "
		<< " BER: " << berc.get_errorrate() << " . "
		<< "Obtained " << per.get_errors() << " error packets. "
		<< " PER: " << per.get_errorrate() << " . "
		<< endl << endl << flush;

	cout << "Done!" << endl << "Success to demodulate the bit stream !" << endl ;
	cout << "Please evaluate the function and performance." << endl ;

	double countIterationAverage = 0.0f;
	for (int i=0;i<COUNT_REPEAT;i++)
	{
		//cout << countIteration[i] << " iteration, " ;

		if (countIteration[i]<0)
			countIteration[i] *= -1;

		countIterationAverage += countIteration[i];
	}
	countIterationAverage = (int)(countIterationAverage/COUNT_REPEAT+0.5) ;
	
	cout << endl << countIterationAverage << " iterations in decoding one ldpc code" << endl << endl ;

	//free( bitOut );
	free( bitsPacketsPadding );


	free( bitsLDPC );
	free( softbits_buf );

	cudaDeviceReset();

	system( "pause" );

	return 0;
}
