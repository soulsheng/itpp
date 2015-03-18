// Generate some example LDPC codes

#include <itpp/itcomm.h>
#include "dvbUtility.h"

using namespace itpp;
using namespace std;

#define		FILENAME_IT		"../data/random_3_6_16200.it"
#define		EBNO			100//2.20

#define		COUNT_REPEAT_DEF	1	// repeat time 
#define		SIZE_PACKET		188

#define		N_BCH			31
#define		T_BCH			2
#define		K_BCH			21

#define		VAR_SIZE_CODE		16200
#define		CHECK_SIZE_CODE		8100//8073

enum	MOD_TYPE
{
	MOD_BPSK,	//	0-
	MOD_QPSK,	//	1-
	MOD_8PSK,	//	2-
	MOD_16APSK,	//	3- 
	MOD_32APSK	//	4-
};

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

	  MOD_TYPE	modType = MOD_QPSK;

#if 0
	  QPSK qpsk;
#else
	  SymbolTable* pSymbol = new SymbolTable(2);
	  Modulator_2D* pModulator = new Modulator_2D( pSymbol->getSymbols(), pSymbol->getBits10Symbols() );
#endif

	  BPSK bpsk;
#if 0
	  APSK32 apsk32(32);
	  APSK16 apsk16(16);
#else
	  SymbolTable* pSymbol4 = new SymbolTable(4);
	  Modulator_2D* pModulator4 = new Modulator_2D( pSymbol->getSymbols(), pSymbol->getBits10Symbols() );

	  SymbolTable* pSymbol5 = new SymbolTable(5);
	  Modulator_2D* pModulator5 = new Modulator_2D( pSymbol->getSymbols(), pSymbol->getBits10Symbols() );
#endif
	  // Noise variance is N0/2 per dimension
	  double N0 = pow(10.0, -EBNO / 10.0) / ldpc.get_rate();
	  AWGN_Channel chan(N0 / 2);


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

	  int nSizeMod = nldpc;

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

		  // step 2: bch encode
		  bvec bitsoutBCHEnc = bch.encode(bitsinBCHEnc);

		  bvec bitsinLDPCEnc = zeros_b(kldpc);
		  bitsinLDPCEnc.set_subvector(0, bitsoutBCHEnc);

		  // step 3: ldpc encode
		  bvec bitsoutLDPCEnc = ldpc.encode(bitsinLDPCEnc);


		  // step 4-6: modulate	-- awgn -- Demodulate
		  vec	dMOD;	// double vector
		  cvec	cMOD;	// complex vector

		  // Received data
		  vec	dAWGN;
		  cvec	cAWGN;

		  switch ( modType )
		  {
		  case MOD_BPSK:
			  dMOD = bpsk.modulate_bits(bitsoutLDPCEnc);
			  dAWGN = chan(dMOD);
			  convertVecToBuffer( bitsMOD, dAWGN );

			  nSizeMod = nldpc*2/bpsk.get_k();

			  break;

		  case MOD_QPSK:

			  cout << "bitsoutLDPCEnc.left(16)" << bitsoutLDPCEnc.left(16) << endl;

			  cMOD = pModulator->modulate_bits(bitsoutLDPCEnc);
			  cout << "cMOD.left(8)" << cMOD.left(8) << endl;

			  cAWGN = chan(cMOD);
			  convertVecToBuffer( bitsMOD, cAWGN );

			  cout << "cMOD.left(8)" << cAWGN.left(8) << endl;

			  nSizeMod = nldpc*2/pModulator->get_k();

			  break;

		  case MOD_16APSK:

			  cout << "bitsoutLDPCEnc.left(16)" << bitsoutLDPCEnc.left(16) << endl;

			  cMOD = pModulator4->modulate_bits(bitsoutLDPCEnc);
			  cout << "cMOD.left(4)" << cMOD.left(4) << endl;
			  cAWGN = chan(cMOD);

			  convertVecToBuffer( bitsMOD, cMOD );

			  cout << "cAWGN.left(4)" << cAWGN.left(4) << endl;

			  nSizeMod = nldpc*2/pModulator4->get_k();

			  break;

		  case MOD_32APSK:

			  cout << "bitsoutLDPCEnc.left(15)" << bitsoutLDPCEnc.left(15) << endl;

			  cMOD = pModulator5->modulate_bits(bitsoutLDPCEnc);
#if 0
			  cAWGN = chan(cMOD);
			  convertVecToBuffer( bitsMOD, cAWGN );
#else
			  convertVecToBuffer( bitsMOD, cMOD );

			  cout << "cMOD.left(3)" << cMOD.left(3) << endl;
#endif

			  nSizeMod = nldpc*2/pModulator5->get_k();

			  break;

		  default:
			  break;;
		  }

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
