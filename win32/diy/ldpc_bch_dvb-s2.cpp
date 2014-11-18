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

int main(int argc, char **argv)
{
	// step 0: intialize ldpc,bch,bpsk,awgn
	LDPC_Generator_Systematic G; // for codes created with ldpc_gen_codes since generator exists
	LDPC_Code ldpc(FILENAME_IT, &G);
	
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

		// step 1: input message
		bvec bitsinBCHEnc = randb(Kbch);

		// step 2: bch encode
		bvec bitsoutBCHEnc(Nbch);
		for (int j = 0; j < nSplit; j++)
		{
			bvec bitsinBCHOne = bitsinBCHEnc.mid(j*K_BCH, K_BCH);
			bvec encodedBCHOne = bch.encode(bitsinBCHOne);
			bitsoutBCHEnc.set_subvector(j*N_BCH, encodedBCHOne);
		}


		bvec bitsinLDPCEnc(kldpc);
		bitsinLDPCEnc.set_subvector(0, bitsoutBCHEnc);

		// step 3: ldpc encode
		bvec bitsoutLDPCEnc = ldpc.encode(bitsinLDPCEnc);
		
		// step 4: modulate
		vec s = Mod.modulate_bits(bitsoutLDPCEnc);

		
		// step 5: awgn Received data
		vec x = chan(s);

		// step 6: Demodulate
		vec softbits = Mod.demodulate_soft_bits(x, N0);

		// step 7: ldpc Decode the received bits
		QLLRvec llr;
		ldpc.bp_decode(ldpc.get_llrcalc().to_qllr(softbits), llr);
		bvec bitsoutLDPCDec = llr < 0;
		bvec bitsinBCHDec = bitsoutLDPCDec.left(Nbch);
		//      bvec bitsout = C.decode(softbits); // (only systematic bits)

		// step 8: bch decode
		bvec bitsoutBCHDec(Kbch);
		for (int j = 0; j < nSplit; j++)
		{
			bvec bitsinBCHDecOne = bitsinBCHDec.mid(j*N_BCH, N_BCH);
			bvec bitsoutBCHDecOne = bch.decode(bitsinBCHDecOne);
			bitsoutBCHDec.set_subvector(j*K_BCH, bitsoutBCHDecOne);
		}

		// step 9: verify result, Count the number of errors
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
