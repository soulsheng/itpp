#include <itpp/itcomm.h>
#include <sstream>

using namespace std;
using namespace itpp;

#define		EBNO			6.3
#define		N_BCH			15
#define		T_BCH			1
#define		K_BCH			11

int main(int argc, char **argv)
{

	BCH C(N_BCH, T_BCH);

	BPSK Mod;

	// Noise variance is N0/2 per dimension
	double N0 = pow(10.0, -EBNO / 10.0) / C.get_rate();
	AWGN_Channel chan(N0 / 2);

	BERC berc;  // Counters for coded and uncoded BER

	RNG_randomize();

	int64_t Nbits = 3960LL; // maximum number of bits simulated 5000000000LL
    for (int64_t i = 0; i < Nbits; i += K_BCH) 
	{
      
		bvec bitsin = randb(K_BCH);

		bvec codbits = C.encode(bitsin);

		vec s = Mod.modulate_bits(codbits);

		
		// Received data
		vec x = chan(s);

		// Demodulate
		bvec decbits = Mod.demodulate_bits(x);

		// Decode the received bits
		bvec bitsout = C.decode(decbits); 

		// Count the number of errors
		berc.count(bitsin, bitsout);

      
		cout << "Eb/N0 = " << EBNO << "  Simulated "
				<< berc.get_total_bits() << " bits. " << endl
				<< "Obtained " << berc.get_errors() << " bit errors. "
				<< " BER: " << berc.get_errorrate() << " . "
				<< endl << flush;
    }


	return 0;
}
