#include <itpp/itcomm.h>
#include <sstream>

using namespace std;
using namespace itpp;

#define		FILENAME_IT		"RU_1000.it"
#define		EBNO			0.67

int main(int argc, char **argv)
{

	LDPC_Generator_Systematic G; // for codes created with ldpc_gen_codes since generator exists
	LDPC_Code C(FILENAME_IT, &G);
	cout << C << endl;

	// High performance: 2500 iterations, high resolution LLR algebra
	C.set_exit_conditions(2500);

	// Alternate high speed settings: 50 iterations, logmax approximation
	// C.set_llrcalc(LLR_calc_unit(12,0,7));


	BPSK Mod;

	// Noise variance is N0/2 per dimension
	double N0 = pow(10.0, -EBNO / 10.0) / C.get_rate();
	AWGN_Channel chan(N0 / 2);

	BERC berc;  // Counters for coded and uncoded BER
	BLERC ferc; // Counter for coded FER
	ferc.set_blocksize(C.get_nvar() - C.get_ncheck());

	RNG_randomize();

	int64_t Nbits = 10000LL; // maximum number of bits simulated 5000000000LL
    for (int64_t i = 0; i < Nbits; i += C.get_nvar()) 
	{
      
		int N = C.get_nvar();             // number of bits per codeword
		bvec bitsin = zeros_b(N);
		vec s = Mod.modulate_bits(bitsin);

		
		// Received data
		vec x = chan(s);

		// Demodulate
		vec softbits = Mod.demodulate_soft_bits(x, N0);

		// Decode the received bits
		QLLRvec llr;
		C.bp_decode(C.get_llrcalc().to_qllr(softbits), llr);
		bvec bitsout = llr < 0;
		//      bvec bitsout = C.decode(softbits); // (only systematic bits)

		// Count the number of errors
		berc.count(bitsin, bitsout);
		ferc.count(bitsin, bitsout);

      
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
