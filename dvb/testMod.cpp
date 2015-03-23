// Generate some example LDPC codes

#include <itpp/itcomm.h>
#include "dvbUtility.h"
#include "modulatorFactory.h"
#include "ldpc_bp_decode.h"
#include "ldpc_bp_decode.cuh"

using namespace itpp;
using namespace std;

#define TEST_HARD_BIT	0
#define	USE_GPU			1

int main(int argc, char **argv)
{
	  // ldpc
	ifstream  testfile;
	testfile.open( FILENAME_IT );
	if ( testfile == NULL )	{
		cout << "Can not find ldpc code file - \"random_3_6_16200.it\" in data path!" << endl << "Please run ldpc_gen_codes.exe to generate one.";
		return 0;
	}	else	{
		cout << "Success to load ldpc code file - \"random_3_6_16200.it\" in data path!" << endl ;
	}
	testfile.close();

	cout << "Generating ldpc data file - \"bitfile.dat\" to data path!" << endl ;
	cout << "Please wait for a few minutes ..." << endl;

	LDPC_Generator_Systematic G; // for codes created with ldpc_gen_codes since generator exists	
	LDPC_Code ldpc(FILENAME_IT, &G);

	  // Noise variance is N0/2 per dimension
	  double N0 = pow(10.0, -EBNO / 10.0) / ldpc.get_rate();
	  AWGN_Channel chan(N0 / 2);

	  ModulatorFactory	mods;	// 调制解调器件库

	  BERC berc;  // Counters for coded and uncoded BER
	  BLERC per(SIZE_PACKET); // Counter for coded FER

	  int nmaxX1 = max(ldpc.sumX1._data(), ldpc.sumX1.size());
	  int nmaxX2 = max(ldpc.sumX2._data(), ldpc.sumX2.size());

	  ldpc_gpu	ldpc_gpu_diy;
	  ldpc_gpu_diy.initialize(ldpc.nvar, ldpc.ncheck, 
		  nmaxX1, nmaxX2, 
		  ldpc.sumX1._data(), ldpc.sumX2._data(), ldpc.iind._data(), ldpc.jind._data(), ldpc.V._data(), 	// Parity check matrix parameterization
		  ldpc.mvc._data(), ldpc.mcv._data(),	// temporary storage for decoder (memory allocated when codec defined)
		  ldpc.llrcalc.Dint1, ldpc.llrcalc.Dint2, ldpc.llrcalc.Dint3,	//! Decoder (lookup-table) parameters
		  ldpc.llrcalc.logexp_table._data());

	  int nldpc = VAR_SIZE_CODE;
	  int kldpc = VAR_SIZE_CODE*ldpc.get_rate();
	  char *bitsPacketsPadding = new char[kldpc];

	  char * llrOut = (char*)malloc( nldpc * sizeof(char) );

	  int COUNT_REPEAT = COUNT_REPEAT_DEF;

	  ivec		countIteration(COUNT_REPEAT);

#if REMOVE_BCH
	  int nLengthMSG = kldpc;
#else
	  int nLengthMSG = Kbch;
#endif

	  for (int64_t i = 0; i < COUNT_REPEAT; i ++) 
	  {
		  // step 0: prepare input packets from rand data or file stream
		  memset( bitsPacketsPadding, 0, sizeof(char)*nLengthMSG );
		  srand( (unsigned int)i*CHECK_SIZE_CODE ) ;

		  int nCountPacket = nLengthMSG / SIZE_PACKET;
		  for (int j = 0; j < nCountPacket; j ++) 
		  {
			  char *onePacket = bitsPacketsPadding + j*SIZE_PACKET;
			  for (int k = 0; k < SIZE_PACKET; k ++) 
			  {
				onePacket[k] = rand()%2;
			  }
		  }

		  // step 1: input message
		  bvec bitsinLDPCEnc( kldpc );
		  convertBufferToVec( bitsPacketsPadding, bitsinLDPCEnc );
		  //cout << " bitsinLDPCEnc.left(16)" << bitsinLDPCEnc.left(16) << endl;

		  // ldpc encode
		  bvec bitsoutLDPCEnc = ldpc.encode(bitsinLDPCEnc);
		  //cout << "bitsoutLDPCEnc.left(16)" << bitsoutLDPCEnc.left(16) << endl;

		  // step 4-6: modulate	-- awgn -- Demodulate

		  int nModTypeRAND = rand()%4 + 2;
		  cout << "nModTypeRAND = " << nModTypeRAND << endl;

		  MOD_TYPE	modType = (MOD_TYPE)nModTypeRAND;
		  Modulator_2D* pModulator = mods.findModulator( modType );

		  if ( NULL == pModulator )
		  {
			  cout << "can not modulate" << endl;
			  continue;
		  }

		  cvec	cMOD = pModulator->modulate_bits(bitsoutLDPCEnc);
			//cout << "cMOD.left(8)" << cMOD.left(8) << endl;

#if REMOVE_NOISE
			cvec	cAWGN = cMOD;
#else
			cvec	cAWGN = chan(cMOD);

			//cout << "cAWGN.left(8)" << cAWGN.left(8) << endl;
#endif

		  int	nSizeMod = nldpc*2/pModulator->get_k();

		  // demodulate
#if TEST_HARD_BIT
		  bvec hardbits = pModulator->demodulate_bits( cAWGN );
		  cout << "hardbits.left(16)" << hardbits.left(16) << endl;
		  vec softbits(nldpc);
#else
		  vec softbits = pModulator->demodulate_soft_bits(cAWGN, N0, APPROX);
		  //cout << "softbits.left(16)" << softbits.left(16) << endl;
#endif

		  QLLRvec llrIn = ldpc.get_llrcalc().to_qllr(softbits);

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
		  bvec bitsoutLDPCDec(nldpc);
		  for (int j=0;j<nldpc;j++)
			  bitsoutLDPCDec[j] = llrOut[j];


		  berc.count(bitsinLDPCEnc, bitsoutLDPCDec);
		  per.count(bitsinLDPCEnc, bitsoutLDPCDec);

		  //cout << " bitsinLDPCEnc.left(16)" << bitsinLDPCEnc.left(16) << endl;
		  //cout << "bitsoutLDPCDec.left(16)" << bitsoutLDPCDec.left(16) << endl;


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

	  double countIterationAverage = 0.0f;
	  for (int i=0;i<COUNT_REPEAT;i++)
	  {
		  cout << countIteration[i] << " iteration, " ;

		  if (countIteration[i]<0)
			  countIteration[i] *= -1;

		  countIterationAverage += countIteration[i];
	  }
	  countIterationAverage = (int)(countIterationAverage/COUNT_REPEAT+0.5) ;
	  cout << endl << countIterationAverage << " iterations in decoding one ldpc code" << endl << endl ;

	  free( bitsPacketsPadding );
	  free( llrOut );

	  system( "pause" );

	return 0;
}
