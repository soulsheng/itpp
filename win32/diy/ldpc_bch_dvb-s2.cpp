#include <itpp/itcomm.h>
#include <sstream>

using namespace std;
using namespace itpp;

#define		FILENAME_IT		"RU_10000.it"
#define		EBNO			1.1
#define		COUNT_REPEAT	1000	// repeat time 

#define		N_BCH			31
#define		T_BCH			2
#define		K_BCH			21

int main(int argc, char **argv)
{

	LDPC_Generator_Systematic G; // for codes created with ldpc_gen_codes since generator exists
	LDPC_Code ldpc(FILENAME_IT, &G);
	cout << ldpc << endl;

	// High performance: 2500 iterations, high resolution LLR algebra
	ldpc.set_exit_conditions(2500);

	// Alternate high speed settings: 50 iterations, logmax approximation
	// C.set_llrcalc(LLR_calc_unit(12,0,7));
	
	BCH bch(N_BCH, T_BCH);


	BPSK Mod;

	// Noise variance is N0/2 per dimension
	double N0 = pow(10.0, -EBNO / 10.0) / ldpc.get_rate();
	AWGN_Channel chan(N0 / 2);

	BERC berc;  // Counters for coded and uncoded BER
	BLERC ferc; // Counter for coded FER

	RNG_randomize();


	int kldpc = ldpc.get_ninfo();             // number of bits per codeword

	int nSplit = kldpc / N_BCH;
	int Nbch = nSplit * N_BCH;

	int Kbch = nSplit * K_BCH;

	ferc.set_blocksize(Kbch);

    for (int64_t i = 0; i < COUNT_REPEAT; i ++) 
	{

		bvec bitsinBCHEnc = randb(Kbch);

		bvec bitsoutBCHEnc(Nbch);
		for (int j = 0; j < nSplit; j++)
		{
			bvec bitsinBCHOne = bitsinBCHEnc.mid(j*K_BCH, K_BCH);
			bvec encodedBCHOne = bch.encode(bitsinBCHOne);
			bitsoutBCHEnc.set_subvector(j*N_BCH, encodedBCHOne);
		}


		bvec bitsinLDPCEnc(kldpc);
		bitsinLDPCEnc.set_subvector(0, bitsoutBCHEnc);

		bvec bitsoutLDPCEnc = ldpc.encode(bitsinLDPCEnc);
		
		vec s = Mod.modulate_bits(bitsoutLDPCEnc);

		
		// Received data
		vec x = chan(s);

		// Demodulate
		vec softbits = Mod.demodulate_soft_bits(x, N0);

		// Decode the received bits
		QLLRvec llr;
		ldpc.bp_decode(ldpc.get_llrcalc().to_qllr(softbits), llr);
		bvec bitsoutLDPCDec = llr < 0;
		bvec bitsinBCHDec = bitsoutLDPCDec.left(Nbch);
		//      bvec bitsout = C.decode(softbits); // (only systematic bits)

		bvec bitsoutBCHDec(Kbch);
		for (int j = 0; j < nSplit; j++)
		{
			bvec bitsinBCHDecOne = bitsinBCHDec.mid(j*N_BCH, N_BCH);
			bvec bitsoutBCHDecOne = bch.decode(bitsinBCHDecOne);
			bitsoutBCHDec.set_subvector(j*K_BCH, bitsoutBCHDecOne);
		}

		// Count the number of errors
		berc.count(bitsinBCHEnc, bitsoutBCHDec);
		ferc.count(bitsinBCHEnc, bitsoutBCHDec);

      
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


	return 0;
}
