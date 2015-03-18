// Generate some example LDPC codes

#include <itpp/itcomm.h>
#include "dvbUtility.h"
#include "modulatorFactory.h"

using namespace itpp;
using namespace std;

#define		FILENAME_IT		"../data/random_3_6_16200.it"
#define		EBNO			2.20

#define		COUNT_REPEAT_DEF	1	// repeat time 
#define		SIZE_PACKET		188

#define		N_BCH			31
#define		T_BCH			2
#define		K_BCH			21

#define		VAR_SIZE_CODE		16200
#define		CHECK_SIZE_CODE		8100//8073

#define		REMOVE_NOISE		0
#define		REMOVE_BCH			1

//! Maximum value of vector
int max(int *v, int N)
{
	int tmp = v[0];
	for (int i = 1; i < N; i++)
		if (v[i] > tmp)
			tmp = v[i];
	return tmp;
}

//! Minimum value of vector
int min(int *v, int N)
{
	int tmp = v[0];
	for (int i = 1; i < N; i++)
		if (v[i] < tmp)
			tmp = v[i];
	return tmp;
}

int main(int argc, char **argv)
{
	// generate input bit and modulated bit
	ifstream  testfile;
	testfile.open( FILENAME_IT );
	if ( testfile == NULL )
	{
		cout << "Can not find ldpc code file - \"random_3_6_16200.it\" in data path!" << endl << "Please run ldpc_gen_codes.exe to generate one.";
		return 0;
	}
	else
	{
		cout << "Success to load ldpc code file - \"random_3_6_16200.it\" in data path!" << endl ;
	}
	testfile.close();

	cout << "Generating ldpc data file - \"bitfile.dat\" to data path!" << endl ;
	cout << "Please wait for a few minutes ..." << endl;

	  LDPC_Generator_Systematic G; // for codes created with ldpc_gen_codes since generator exists
	  LDPC_Code ldpc(FILENAME_IT, &G);

	  BCH bch(N_BCH, T_BCH);

	  // Noise variance is N0/2 per dimension
	  double N0 = pow(10.0, -EBNO / 10.0) / ldpc.get_rate();
	  AWGN_Channel chan(N0 / 2);

	  ModulatorFactory	mods;	// 调制解调器件库

	  ofstream  bitfile;
	  bitfile.open( "../data/bitfile.dat" );
	  if ( bitfile == NULL )
	  {
		  return 0;
	  }

	  ofstream  bitfileMOD;
	  bitfileMOD.open( "../data/bitfileMOD.dat" );
	  if ( bitfileMOD == NULL )
	  {
		  return 0;
	  }

	  int nldpc = VAR_SIZE_CODE;
	  int kldpc = CHECK_SIZE_CODE;             // number of bits per codeword

	  int nSplit = kldpc / N_BCH;
	  int Nbch = nSplit * N_BCH;

	  int Kbch = nSplit * K_BCH;

	  char *bitsPacketsPadding = new char[Kbch];
	  char *bitsBCH = new char[kldpc];
	  char *bitsLDPC = new char[nldpc];
	  double *bitsMOD = new double[nldpc];

	  int COUNT_REPEAT = COUNT_REPEAT_DEF;
	  bitfile.write( (char*)&COUNT_REPEAT, sizeof(int)*1);
	  bitfile.write( (char*)&Kbch, sizeof(int)*1);

	  for (int64_t i = 0; i < COUNT_REPEAT; i ++) 
	  {
		  // step 0: prepare input packets from rand data or file stream
		  memset( bitsPacketsPadding, 0, sizeof(char)*Kbch );
		  srand( (unsigned int)i*CHECK_SIZE_CODE ) ;

		  int nCountPacket = Kbch / SIZE_PACKET;
		  for (int j = 0; j < nCountPacket; j ++) 
		  {
			  char *onePacket = bitsPacketsPadding + j*SIZE_PACKET;
			  for (int k = 0; k < SIZE_PACKET; k ++) 
			  {
				onePacket[k] = rand()%2;
			  }
		  }

		  bitfile.write(bitsPacketsPadding, sizeof(char)*Kbch);


		  // step 1: input message
		  bvec bitsinBCHEnc( Kbch );
		  convertBufferToVec( bitsPacketsPadding, bitsinBCHEnc );

		  bvec bitsinLDPCEnc = zeros_b(kldpc);

#if	REMOVE_BCH
		  bitsinLDPCEnc.set_subvector(0, bitsinBCHEnc);
#else
		  // step 2: bch encode
		  bvec bitsoutBCHEnc = bch.encode(bitsinBCHEnc);

		  bitsinLDPCEnc.set_subvector(0, bitsoutBCHEnc);
#endif

		  // step 3: ldpc encode
		  bvec bitsoutLDPCEnc = ldpc.encode(bitsinLDPCEnc);
		  cout << "bitsoutLDPCEnc.left(16)" << bitsoutLDPCEnc.left(16) << endl;

		  // step 4-6: modulate	-- awgn -- Demodulate
		  MOD_TYPE	modType = MOD_32APSK;
		  Modulator_2D* pModulator = mods.findModulator( modType );

		  if ( NULL == pModulator )
		  {
			  cout << "can not modulate" << endl;
			  continue;
		  }

		  cvec	cMOD = pModulator->modulate_bits(bitsoutLDPCEnc);
			cout << "cMOD.left(8)" << cMOD.left(8) << endl;

#if REMOVE_NOISE
			convertVecToBuffer( bitsMOD, cMOD );	
#else
			cvec	cAWGN = chan(cMOD);
			convertVecToBuffer( bitsMOD, cAWGN );

			cout << "cAWGN.left(8)" << cAWGN.left(8) << endl;
#endif

		  int	nSizeMod = nldpc*2/pModulator->get_k();

			
		  for ( int j=0; j<nSizeMod; j++ )
		  	  bitfileMOD << bitsMOD[j] << " " ;

	  }

	  bitfile.close();
	  bitfileMOD.close();

	  free( bitsPacketsPadding );
	  free( bitsMOD );

	  if ( bitfile != NULL )
	  {
		  cout << "Done!" << endl << "Success to generate ldpc data file - \"bitfile.dat\" in data path!" << endl ;
		  cout << "Please run dvb-s2.exe to demodulate and decode it." << endl ;
		  return 0;
	  }

	  system( "pause" );

	return 0;
}
