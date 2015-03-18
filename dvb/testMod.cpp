// Generate some example LDPC codes

#include <itpp/itcomm.h>
#include "dvbUtility.h"
#include "modulatorFactory.h"

using namespace itpp;
using namespace std;

#define TEST_HARD_BIT	1

int main(int argc, char **argv)
{

	  // Noise variance is N0/2 per dimension
	  double N0 = pow(10.0, -EBNO / 10.0) / 0.5f;
	  AWGN_Channel chan(N0 / 2);

	  ModulatorFactory	mods;	// 调制解调器件库

	  BERC berc;  // Counters for coded and uncoded BER
	  BLERC per(SIZE_PACKET); // Counter for coded FER

	  int nldpc = VAR_SIZE_CODE;

	  char *bitsPacketsPadding = new char[nldpc];

	  int COUNT_REPEAT = COUNT_REPEAT_DEF;

	  for (int64_t i = 0; i < COUNT_REPEAT; i ++) 
	  {
		  // step 0: prepare input packets from rand data or file stream
		  memset( bitsPacketsPadding, 0, sizeof(char)*nldpc );
		  srand( (unsigned int)i*CHECK_SIZE_CODE ) ;

		  int nCountPacket = nldpc / SIZE_PACKET;
		  for (int j = 0; j < nCountPacket; j ++) 
		  {
			  char *onePacket = bitsPacketsPadding + j*SIZE_PACKET;
			  for (int k = 0; k < SIZE_PACKET; k ++) 
			  {
				onePacket[k] = rand()%2;
			  }
		  }

		  // step 1: input message
		  bvec bitsin( nldpc );
		  convertBufferToVec( bitsPacketsPadding, bitsin );
		  cout << "bitsin.left(16)" << bitsin.left(16) << endl;


		  // step 4-6: modulate	-- awgn -- Demodulate
		  MOD_TYPE	modType = (MOD_TYPE)MOD_TYPE_DEFAULT;
		  Modulator_2D* pModulator = mods.findModulator( modType );

		  if ( NULL == pModulator )
		  {
			  cout << "can not modulate" << endl;
			  continue;
		  }

		  cvec	cMOD = pModulator->modulate_bits(bitsin);
			cout << "cMOD.left(8)" << cMOD.left(8) << endl;

#if REMOVE_NOISE
			cvec	cAWGN = cMOD;
#else
			cvec	cAWGN = chan(cMOD);

			cout << "cAWGN.left(8)" << cAWGN.left(8) << endl;
#endif

		  int	nSizeMod = nldpc*2/pModulator->get_k();

		  // demodulate
#if TEST_HARD_BIT
		  bvec hardbits = pModulator->demodulate_bits( cAWGN );
		  cout << "hardbits.left(16)" << hardbits.left(16) << endl;
		  vec softbits(nldpc);
#else
		  vec softbits = pModulator->demodulate_soft_bits(cAWGN, N0);
		  cout << "softbits.left(16)" << softbits.left(16) << endl;
#endif


		  berc.count(bitsin, hardbits);
		  per.count(bitsin, hardbits);

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

	  free( bitsPacketsPadding );

	  system( "pause" );

	return 0;
}
