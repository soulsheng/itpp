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

enum	STEP_TIMER
{
	STEP_ENC_BCH,	//	0-
	STEP_ENC_LDPC,	//	1-
	STEP_MOD,		//	2- MOD AWGN DEMOD
	STEP_DEC_LDPC,	//	3-
	STEP_DEC_BCH,	//	4-
	COUNT_TIMER_STEP
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

	Real_Timer	timer, timerStep[COUNT_TIMER_STEP];
	vec			timerValue(COUNT_REPEAT);

	Vec<vec>			timerStepValue(COUNT_TIMER_STEP);
	for (int i=0;i<COUNT_TIMER_STEP;i++)
		timerStepValue[i].set_size(COUNT_REPEAT);


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

		int nStep = STEP_ENC_BCH;
		timerStep[nStep].reset();
		timerStep[nStep].start();

		// step 2: bch encode
		bvec bitsoutBCHEnc = zeros_b(Nbch);
		for (int j = 0; j < nSplit; j++)
		{
			bvec bitsinBCHOne = bitsinBCHEnc.mid(j*K_BCH, K_BCH);
			bvec encodedBCHOne = bch.encode(bitsinBCHOne);
			bitsoutBCHEnc.set_subvector(j*N_BCH, encodedBCHOne);
		}

		timerStep[nStep].stop();
		timerStepValue[nStep][i] = timerStep[nStep].get_time() ;

		bvec bitsinLDPCEnc = zeros_b(kldpc);
		bitsinLDPCEnc.set_subvector(0, bitsoutBCHEnc);

		nStep ++ ;
		timerStep[nStep].reset();
		timerStep[nStep].start();

		// step 3: ldpc encode
		bvec bitsoutLDPCEnc = ldpc.encode(bitsinLDPCEnc);
		
		timerStep[nStep].stop();
		timerStepValue[nStep][i] = timerStep[nStep].get_time() ;


		nStep ++ ;
		timerStep[nStep].reset();
		timerStep[nStep].start();

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


		timerStep[nStep].stop();
		timerStepValue[nStep][i] = timerStep[nStep].get_time() ;


		nStep ++ ;
		timerStep[nStep].reset();
		timerStep[nStep].start();

		// step 7: ldpc Decode the received bits
		QLLRvec llr;
		ldpc.bp_decode(ldpc.get_llrcalc().to_qllr(softbits), llr);
		bvec bitsoutLDPCDec = llr < 0;
		bvec bitsinBCHDec = bitsoutLDPCDec.left(Nbch);
		//      bvec bitsout = C.decode(softbits); // (only systematic bits)


		timerStep[nStep].stop();
		timerStepValue[nStep][i] = timerStep[nStep].get_time() ;

		nStep ++ ;
		timerStep[nStep].reset();
		timerStep[nStep].start();

		// step 8: bch decode
		bvec bitsoutBCHDec = zeros_b(Kbch);
		for (int j = 0; j < nSplit; j++)
		{
			bvec bitsinBCHDecOne = bitsinBCHDec.mid(j*N_BCH, N_BCH);
			bvec bitsoutBCHDecOne = bch.decode(bitsinBCHDecOne);
			bitsoutBCHDec.set_subvector(j*K_BCH, bitsoutBCHDecOne);
		}

		timerStep[nStep].stop();
		timerStepValue[nStep][i] = timerStep[nStep].get_time() ;

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

	double timerAverageAll = 0.0f, timerStepAverage[COUNT_TIMER_STEP]={0};
	for (int i=0;i<COUNT_REPEAT;i++)
	{
		cout << timerValue[i] << " s, " ;
		timerAverageAll += timerValue[i];

		for (int j=0;j<COUNT_TIMER_STEP;j++)
			timerStepAverage[j] += timerStepValue[j][i];
	}

	cout << endl << timerAverageAll/COUNT_REPEAT << " s Average all" << endl ;
	
	for (int j=0;j<COUNT_TIMER_STEP;j++)
		cout << endl << timerStepAverage[j]/COUNT_REPEAT << " s Average step " << STEP_TIMER(j) << endl ;

	return 0;
}
