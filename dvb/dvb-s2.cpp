#include <itpp/itcomm.h>
#include <sstream>
#include <cuda_runtime.h>
#include "ldpc_bp_decode.h"
#include "ldpc_bp_decode.cuh"
#include "helper_timer.h"
#include "driverUtility.h"

using namespace std;
using namespace itpp;

//#define		FILENAME_IT		"../data/RU_16200.it"
#define		FILENAME_IT		"../data/random_3_6_16200.it"
#define		EBNO			2.15
#define		COUNT_REPEAT	100	// repeat time 

#define		N_BCH			31
#define		T_BCH			2
#define		K_BCH			21

#define		SIZE_PACKET		188

#define		USE_GPU		1

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
	// step 0: intialize ldpc,bch,bpsk,awgn
	LDPC_Generator_Systematic G; // for codes created with ldpc_gen_codes since generator exists
	LDPC_Code ldpc(FILENAME_IT, &G);
	
	BCH bch(N_BCH, T_BCH);

	MOD_TYPE	modType = MOD_QPSK;

	QPSK qpsk;
	BPSK bpsk;

	// Noise variance is N0/2 per dimension
	double N0 = pow(10.0, -EBNO / 10.0) / ldpc.get_rate();
	AWGN_Channel chan(N0 / 2);

	BERC berc;  // Counters for coded and uncoded BER
	BLERC per(SIZE_PACKET); // Counter for coded FER

	RNG_randomize();

	StopWatchInterface	*timer, *timerStep;
	sdkCreateTimer( &timer );
	sdkCreateTimer( &timerStep );

	vec			timerValue(COUNT_REPEAT);

	vec			timerStepValue(COUNT_REPEAT);
	ivec		countIteration(COUNT_REPEAT);

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

    for (int64_t i = 0; i < COUNT_REPEAT; i ++) 
	{
		// step 0: prepare input packets from rand data or file stream
		bvec	bitsPacketsPadding = zeros_b(Kbch);
		int nCountPacket = Kbch / SIZE_PACKET;
		for (int j = 0; j < nCountPacket; j ++) 
		{
			bvec onePacket = randb(SIZE_PACKET);
			bitsPacketsPadding.set_subvector(j*SIZE_PACKET, onePacket);
		}

		sdkResetTimer( &timer );
		sdkStartTimer( &timer );

		// step 1: input message
		bvec bitsinBCHEnc = bitsPacketsPadding;//randb(Kbch);

		// step 2: bch encode
		bvec bitsoutBCHEnc = zeros_b(Nbch);
		for (int j = 0; j < nSplit; j++)
		{
			bvec bitsinBCHOne = bitsinBCHEnc.mid(j*K_BCH, K_BCH);
			bvec encodedBCHOne = bch.encode(bitsinBCHOne);
			bitsoutBCHEnc.set_subvector(j*N_BCH, encodedBCHOne);
		}

		bvec bitsinLDPCEnc = zeros_b(kldpc);
		bitsinLDPCEnc.set_subvector(0, bitsoutBCHEnc);


		// step 3: ldpc encode
		bvec bitsoutLDPCEnc = ldpc.encode(bitsinLDPCEnc);
		
		// step 4-6: modulate	-- awgn -- Demodulate
		vec		dMOD;	// double vector
		cvec	cMOD;	// complex vector

		// Received data
		vec		dAWGN;
		cvec	cAWGN;

		// Demodulate
		vec softbits;
	
		switch ( modType )
		{
		case MOD_BPSK:
			dMOD = bpsk.modulate_bits(bitsoutLDPCEnc);
			dAWGN = chan(dMOD);
			softbits = bpsk.demodulate_soft_bits(dAWGN, N0);
			break;

		case MOD_QPSK:
			cMOD = qpsk.modulate_bits(bitsoutLDPCEnc);
			cAWGN = chan(cMOD);
			softbits = qpsk.demodulate_soft_bits(cAWGN, N0);
			break;

		default:
			break;;
		}




		// step 7: ldpc Decode the received bits
		//QLLRvec llr(nldpc);
		QLLRvec llrIn = ldpc.get_llrcalc().to_qllr(softbits);
		
#if WRITE_FILE_FOR_DRIVER
		if( i==0 )
			writeArray( llrIn._data(), ldpc.nvar, "../data/input.txt" );
#endif

		sdkResetTimer( &timerStep );
		sdkStartTimer( &timerStep );

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
		timerStepValue[i] = sdkGetTimerValue( &timerStep );


		bvec bitsoutLDPCDec(nldpc);
		for (int j=0;j<nldpc;j++)
			bitsoutLDPCDec[j] = llrOut[j];

		bvec bitsinBCHDec = bitsoutLDPCDec.left(Nbch);
		//      bvec bitsout = C.decode(softbits); // (only systematic bits)


		// step 8: bch decode
		bvec bitsoutBCHDec = zeros_b(Kbch);
		for (int j = 0; j < nSplit; j++)
		{
			bvec bitsinBCHDecOne = bitsinBCHDec.mid(j*N_BCH, N_BCH);
			bvec bitsoutBCHDecOne = bch.decode(bitsinBCHDecOne);
			bitsoutBCHDec.set_subvector(j*K_BCH, bitsoutBCHDecOne);
		}

		// step 9: verify result, Count the number of errors
		berc.count(bitsinBCHEnc, bitsoutBCHDec);
		per.count(bitsinBCHEnc, bitsoutBCHDec);

		sdkStopTimer( &timer );
		timerValue[i] = sdkGetTimerValue( &timer );
        
		cout << "Eb/N0 = " << EBNO << "  Simulated "
				<< i << " frames and "
				<< berc.get_total_bits() << " bits. " << endl
				<< "Obtained " << berc.get_errors() << " bit errors. "
				<< " BER: " << berc.get_errorrate() << " . "
				<< per.get_total_blocks() << " packets. " << endl
				<< "Obtained " << per.get_errors() << " error packets. "
				<< " PER: " << per.get_errorrate() << " . "
				<< endl << flush;
    }

	double timerAverageAll = 0.0f, timerStepAverage = 0.0f;
	for (int i=0;i<COUNT_REPEAT;i++)
	{
		cout << timerValue[i] << " ms, " ;
		timerAverageAll += timerValue[i];
	}
	cout << endl << endl ;

	for (int i=0;i<COUNT_REPEAT;i++)
	{
		cout  << "timerStepValue[ " << i << " ] = "<< timerStepValue[i] << " ms, " << endl;
		timerStepAverage += timerStepValue[i];
	}
	cout << endl << endl ;

	double countIterationAverage = 0.0f;
	for (int i=0;i<COUNT_REPEAT;i++)
	{
		cout << countIteration[i] << " iteration, " ;

		if (countIteration[i]<0)
			countIteration[i] *= -1;

		countIterationAverage += countIteration[i];
	}
	cout << endl << endl ;

	cout << endl << "DVB S-2 totally costs time: "<< timerAverageAll/COUNT_REPEAT << " ms for each code with length of 16200" << endl ;
	
	cout << endl << timerStepAverage/COUNT_REPEAT << " ms for decoding ldpc" << endl ;
	
	cout << endl << (int)(countIterationAverage/COUNT_REPEAT+0.5) << " iterations in decode ldpc" << endl << endl ;

	free( llrOut );
	sdkDeleteTimer( &timer );
	sdkDeleteTimer( &timerStep );

	cudaDeviceReset();

	return 0;
}
