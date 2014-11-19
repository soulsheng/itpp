#include <itpp/itcomm.h>
#include <sstream>

using namespace std;
using namespace itpp;

#define		FILENAME_IT		"RU_16200.it"
#define		EBNO			0.8
#define		COUNT_REPEAT	100	// repeat time 

#define		N_BCH			31
#define		T_BCH			2
#define		K_BCH			21

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
	BLERC ferc; // Counter for coded FER

	RNG_randomize();

	Real_Timer	timer, timerStep;
	vec			timerValue(COUNT_REPEAT);

	vec			timerStepValue(COUNT_REPEAT);


	int kldpc = ldpc.get_ninfo();             // number of bits per codeword

	int nSplit = kldpc / N_BCH;
	int Nbch = nSplit * N_BCH;

	int Kbch = nSplit * K_BCH;

	ferc.set_blocksize(Kbch);

    for (int64_t i = 0; i < COUNT_REPEAT; i ++) 
	{
		timer.reset();
		timer.start();

		// step 1: input message
		bvec bitsinBCHEnc = randb(Kbch);

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



		timerStep.reset();
		timerStep.start();

		// step 7: ldpc Decode the received bits
		QLLRvec llr;
		ldpc.bp_decode(ldpc.get_llrcalc().to_qllr(softbits), llr);
		bvec bitsoutLDPCDec = llr < 0;
		bvec bitsinBCHDec = bitsoutLDPCDec.left(Nbch);
		//      bvec bitsout = C.decode(softbits); // (only systematic bits)


		timerStep.stop();
		timerStepValue[i] = timerStep.get_time() ;

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
		ferc.count(bitsinBCHEnc, bitsoutBCHDec);

		timer.stop();
		timerValue[i] = timer.get_time() ;
        
		cout << "Eb/N0 = " << EBNO << "  Simulated "
				<< ferc.get_total_blocks() << " frames and "
				<< berc.get_total_bits() << " bits. " << endl
				<< "Obtained " << berc.get_errors() << " bit errors. "
				<< " BER: " << berc.get_errorrate() << " . "
				<< " FER: " << ferc.get_errorrate() << " . "
				<< ferc.get_errors() << " bit errors. "
				<< ferc.get_corrects() << " bit corrects. "
				<< endl << flush;
    }

	double timerAverageAll = 0.0f, timerStepAverage = 0.0f;
	for (int i=0;i<COUNT_REPEAT;i++)
	{
		cout << timerValue[i] << " s, " ;
		timerAverageAll += timerValue[i];

		timerStepAverage += timerStepValue[i];
	}

	cout << endl << timerAverageAll/COUNT_REPEAT << " s Average all" << endl ;
	
	cout << endl << timerStepAverage/COUNT_REPEAT << " s Average step decode ldpc" << endl ;

	return 0;
}
