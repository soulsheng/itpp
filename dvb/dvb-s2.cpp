#include <itpp/itcomm.h>
#include <sstream>
#include <cuda_runtime.h>
#include "ldpc_bp_decode.h"
#include "ldpc_bp_decode.cuh"
#include "helper_timer.h"
//#include "driverUtility.h"
#include "dvbUtility.h"
#include "modAPSK32.h"
#include "modAPSK16.h"

using namespace std;
using namespace itpp;

//#define		FILENAME_IT		"../data/RU_16200.it"
#define		FILENAME_IT		"../data/random_3_6_16200.it"
#define		FILENAME_ALIST	"../data/dvbs2_r12.alist"
#define		EBNO			10//2.20
//#define		COUNT_REPEAT	10	// repeat time 
#define		TIME_STEP		4	

#define		N_BCH			31
#define		T_BCH			2
#define		K_BCH			21

#define		SIZE_PACKET		188

#define		USE_GPU		1
#define		USE_ALIST		0

enum	MOD_TYPE
{
	MOD_BPSK,	//	0-
	MOD_QPSK,	//	1-
	MOD_8PSK,	//	2-
	MOD_16APSK,	//	3- 
	MOD_32APSK	//	4-
};

int main(int argc, char **argv)
{
	ifstream  testfile;
#if USE_ALIST
	testfile.open( FILENAME_ALIST );
#else
	testfile.open( FILENAME_IT );
#endif
	if ( testfile == NULL )
	{
		cout << "Can not find ldpc code file - \"random_3_6_16200.it\" in data path!" << endl ;
		cout << "Please run ldpc_gen_codes.exe to generate one." << endl;
		return 0;
	}
	else
	{
		cout << "Success to load ldpc code file - alist in data path!" << endl ;
	}
	testfile.close();

	// step 0: intialize ldpc,bch,bpsk,awgn
#if USE_ALIST
	LDPC_Parity	lp( FILENAME_ALIST, "alist" );
	LDPC_Code	ldpc( &lp );
#else
	LDPC_Generator_Systematic G; // for codes created with ldpc_gen_codes since generator exists
	LDPC_Code ldpc(FILENAME_IT, &G);
#endif
	BCH bch(N_BCH, T_BCH);

	MOD_TYPE	modType = MOD_QPSK;

#if 0
	QPSK qpsk;
#else
	SymbolTable* pSymbol = new SymbolTable(2);
	Modulator_2D * pModulator = new Modulator_2D( pSymbol->getSymbols(), pSymbol->getBits10Symbols() );
#endif

	BPSK bpsk;
	APSK32 apsk32(32);
	APSK16 apsk16(16);

	// Noise variance is N0/2 per dimension
	double N0 = pow(10.0, -EBNO / 10.0) / ldpc.get_rate();
	AWGN_Channel chan(N0 / 2);

	BERC berc;  // Counters for coded and uncoded BER
	BLERC per(SIZE_PACKET); // Counter for coded FER

	RNG_randomize();

	StopWatchInterface	*timer, *timerStep;
	sdkCreateTimer( &timer );
	sdkCreateTimer( &timerStep );

	int nldpc = ldpc.get_nvar();
	int kldpc = ldpc.get_ninfo();             // number of bits per codeword

	int nSplit = kldpc / N_BCH;
	int Nbch = nSplit * N_BCH;

	int Kbch = nSplit * K_BCH;

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

	ldpc_gpu	ldpc_gpu_diy;
	ldpc_gpu_diy.initialize(ldpc.nvar, ldpc.ncheck, 
		nmaxX1, nmaxX2, 
		ldpc.sumX1._data(), ldpc.sumX2._data(), ldpc.iind._data(), ldpc.jind._data(), ldpc.V._data(), 	// Parity check matrix parameterization
		ldpc.mvc._data(), ldpc.mcv._data(),	// temporary storage for decoder (memory allocated when codec defined)
		ldpc.llrcalc.Dint1, ldpc.llrcalc.Dint2, ldpc.llrcalc.Dint3,	//! Decoder (lookup-table) parameters
		ldpc.llrcalc.logexp_table._data());

#if WRITE_FILE_FOR_DRIVER
	writeArray( ldpc.sumX1._data(), ldpc.nvar, "../data/sumX1.txt" );
	writeArray( ldpc.sumX2._data(), ldpc.ncheck, "../data/sumX2.txt" );

	writeArray( ldpc.iind._data(), ldpc.nvar * nmaxX1, "../data/iind.txt" );
	writeArray( ldpc.jind._data(), ldpc.ncheck * nmaxX2, "../data/jind.txt" );
	writeArray( ldpc.llrcalc.logexp_table._data(), ldpc.llrcalc.Dint2, "../data/logexp.txt" );
#endif

	char * llrOut = (char*)malloc( nldpc * sizeof(char) );

	ifstream  bitfile;
	bitfile.open( "../data/bitfile.dat" );
	if ( bitfile == NULL )
	{
		cout << "Can not find ldpc data file - \"bitfile.dat\" in data path!" << endl ;
		cout << "Please run ldpc_gen_datas.exe to generate one." << endl;
		return 0;
	}
	else
	{
		cout << "Success to load ldpc data file - \"bitfile.dat\" in data path!" << endl << endl;
	}

	ifstream  bitfileMOD;
	bitfileMOD.open( "../data/bitfileMOD.dat" );
	if ( bitfileMOD == NULL )
	{
		return 0;
	}

	int COUNT_REPEAT = 10;
	bitfile.read( (char*)&COUNT_REPEAT, sizeof(int)*1);
	bitfile.read( (char*)&Kbch, sizeof(int)*1);
	//cout << "COUNT_REPEAT = " << COUNT_REPEAT << endl;	// COUNT_REPEAT = 100

	int nSizeMod = nldpc;
	switch( modType )
	{
	case MOD_QPSK:
		nSizeMod = nldpc*2/pModulator->get_k();
		break;

	case MOD_8PSK:
		nSizeMod = nldpc*2/3;
		break;

	case MOD_16APSK:
		nSizeMod = nldpc*2/apsk16.get_k();
		break;

	case MOD_32APSK:
		nSizeMod = nldpc*2/apsk32.get_k();
		break;


	}
	

	char *bitsPacketsPadding = new char[Kbch];
	double *bitsMOD_N  = new double[nSizeMod*COUNT_REPEAT];
	char *bitsLDPC = new char[nldpc];
	double *softbits_buf = new double[nldpc];

	vec			timerValue(COUNT_REPEAT);

	vec			timerStepValue(COUNT_REPEAT*TIME_STEP);
	ivec		countIteration(COUNT_REPEAT);

	bvec bitsinBCHEnc( Kbch );
	bvec bitsinBCHEnc_N( Kbch*COUNT_REPEAT );


	cout << "Demodulating and decoding the bit stream !" << endl ;
	cout << "Please wait for a few seconds ..." << endl << endl;

	for (int64_t i = 0; i < COUNT_REPEAT; i ++) 
	{
		// step 0: prepare input packets from rand data or file stream
		bitfile.read(bitsPacketsPadding, sizeof(char)*Kbch);
		convertBufferToVec( bitsPacketsPadding, bitsinBCHEnc );
		bitsinBCHEnc_N.set_subvector(i*Kbch, bitsinBCHEnc);

		
		// Received data stream, read from file by trunk length nldpc=16200
		for ( int j=0; j<nSizeMod; j++ )
			bitfileMOD >> bitsMOD_N[j+i*nSizeMod] ;
	}

	int nTimeStep = 0;
	for (int64_t i = 0; i < COUNT_REPEAT; i ++) 
	{

		sdkResetTimer( &timer );
		sdkStartTimer( &timer );


		sdkResetTimer( &timerStep );
		sdkStartTimer( &timerStep );

		// demodulate
		vec		dAWGN( nldpc );
		cvec	cAWGN( nSizeMod/2 );

		vec softbits;
		bvec hardbits;
	
		switch ( modType )
		{
		case MOD_BPSK:
			convertBufferToVec( bitsMOD_N+i*nSizeMod, dAWGN );
			softbits = bpsk.demodulate_soft_bits(dAWGN, N0);
			break;

		case MOD_QPSK:
			convertBufferToVec( bitsMOD_N+i*nSizeMod, cAWGN );

			cout << "cAWGN.left(8)" << cAWGN.left(8) << endl;

			softbits = pModulator->demodulate_soft_bits(cAWGN, N0);

			cout << "softbits.left(16)" << softbits.left(16) << endl;

			break;

		case MOD_16APSK:
			convertBufferToVec( bitsMOD_N+i*nSizeMod, cAWGN );

			cout << "cAWGN.left(4)" << cAWGN.left(4) << endl;
#if 0
			hardbits = apsk16.demodulate_bits( cAWGN );
			cout << "hardbits.left(16)" << hardbits.left(16) << endl;
#else
			softbits = apsk16.demodulate_soft_bits(cAWGN, N0);
			cout << "softbits.left(16)" << softbits.left(16) << endl;
#endif

			break;

		case MOD_32APSK:
			convertBufferToVec( bitsMOD_N+i*nSizeMod, cAWGN );
			cout << "cAWGN.left(3)" << cAWGN.left(3) << endl;

#if 1
			hardbits = apsk32.demodulate_bits( cAWGN );
#else
			softbits = pModulator->demodulate_soft_bits(cAWGN, N0);
#endif
			cout << "hardbits.left(15)" << hardbits.left(15) << endl;

			break;

		default:
			break;;
		}


		sdkStopTimer( &timerStep );
		timerStepValue[nTimeStep++] = sdkGetTimerValue( &timerStep );


		sdkResetTimer( &timerStep );
		sdkStartTimer( &timerStep );

		// ldpc Decode the received bits
		//QLLRvec llr(nldpc);
		QLLRvec llrIn = ldpc.get_llrcalc().to_qllr(softbits);
		
#if WRITE_FILE_FOR_DRIVER
		if( i==0 )
			writeArray( llrIn._data(), ldpc.nvar, "../data/input.txt" );
#endif

#if		USE_GPU
		countIteration[i] = ldpc_gpu_diy.bp_decode_once( llrIn._data(), llrOut ); 
#else
		countIteration[i] = bp_decode( llrIn._data(), llrOut, 
			ldpc.nvar, ldpc.ncheck, 
			nmaxX1, nmaxX2, 
			ldpc.V._data(), ldpc.sumX1._data(), ldpc.sumX2._data(), ldpc.iind._data(), ldpc.jind._data(),	// Parity check matrix parameterization
			ldpc.mvc._data(), ldpc.mcv._data(),	// temporary storage for decoder (memory allocated when codec defined)
			//ldpc.llrcalc );		//!< LLR calculation unit
			ldpc.llrcalc.Dint1, ldpc.llrcalc.Dint2, ldpc.llrcalc.Dint3,	//! Decoder (lookup-table) parameters
			ldpc.llrcalc.logexp_table._data());		//! The lookup tables for the decoder

#endif

		sdkStopTimer( &timerStep );
		timerStepValue[nTimeStep++] = sdkGetTimerValue( &timerStep );


		sdkResetTimer( &timerStep );
		sdkStartTimer( &timerStep );

		bvec bitsoutLDPCDec(nldpc);
		for (int j=0;j<nldpc;j++)
			bitsoutLDPCDec[j] = llrOut[j];

		bvec bitsinBCHDec = bitsoutLDPCDec.left(Nbch);


		// step 8: bch decode
		bvec bitsoutBCHDec = bch.decode(bitsinBCHDec);

		sdkStopTimer( &timerStep );
		timerStepValue[nTimeStep++] = sdkGetTimerValue( &timerStep );

		sdkResetTimer( &timerStep );
		sdkStartTimer( &timerStep );

		// step 9: verify result, Count the number of errors
		bitsinBCHEnc = bitsinBCHEnc_N.mid(i*Kbch, Kbch);

		berc.count(bitsinBCHEnc, bitsoutBCHDec);
		per.count(bitsinBCHEnc, bitsoutBCHDec);

		sdkStopTimer( &timerStep );
		timerStepValue[nTimeStep++] = sdkGetTimerValue( &timerStep );

		sdkStopTimer( &timer );
		timerValue[i] = sdkGetTimerValue( &timer );
        
		
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

	double timerAverageAll = 0.0f, timerStepAverage = 0.0f;
	for (int i=0;i<COUNT_REPEAT;i++)
	{
		//cout << timerValue[i] << " ms, " ;
		timerAverageAll += timerValue[i];
	}
	timerAverageAll /= COUNT_REPEAT;
	cout << endl << endl ;

	for (int i=0;i<COUNT_REPEAT*TIME_STEP;i++)
	{
		//cout  << "timerStepValue[ " << i << " ] = "<< timerStepValue[i] << " ms, " << endl;
	}

	for (int i=0;i<COUNT_REPEAT;i++)
	{
		timerStepAverage += timerStepValue[i*TIME_STEP+1];
	}
	timerStepAverage /= COUNT_REPEAT;

	double countIterationAverage = 0.0f;
	for (int i=0;i<COUNT_REPEAT;i++)
	{
		//cout << countIteration[i] << " iteration, " ;

		if (countIteration[i]<0)
			countIteration[i] *= -1;

		countIterationAverage += countIteration[i];
	}
	countIterationAverage = (int)(countIterationAverage/COUNT_REPEAT+0.5) ;

	cout << endl << "DVB S-2 totally costs time: "<< timerAverageAll << " ms for each code with length of 16200" << endl ;
	
	cout << endl << timerStepAverage << " ms for decoding one ldpc code" << endl ;
	
	cout << endl << countIterationAverage << " iterations in decoding one ldpc code" << endl << endl ;

	free( llrOut );
	free( bitsPacketsPadding );
	bitfile.close();

	free( bitsMOD_N );
	bitfileMOD.close();

	free( bitsLDPC );
	free( softbits_buf );

	sdkDeleteTimer( &timer );
	sdkDeleteTimer( &timerStep );

	cudaDeviceReset();

	system( "pause" );

	return 0;
}
