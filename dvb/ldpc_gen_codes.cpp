// Generate some example LDPC codes

#include <itpp/itcomm.h>

using namespace itpp;
using namespace std;

#define		FILENAME_IT		"../data/random_3_6_16200.it"

int main(int argc, char **argv)
{

	cout << "Generating ldpc code file - \"random_3_6_16200.it\" to data path!" << endl ;
	cout << "Please wait for a few hours(may be 3h) ..." << endl << endl;

  if( 1 )
  { // 16200 bits REGULAR (takes 3 hours to run), (3,6) code with k = 8100
    cout << "========= RANDOM (3,6) CODE ==========" << endl;
    LDPC_Parity_Regular H;
    H.generate(16200, 3, 6,
               "rand",  // random unstructured matrix
               "500 10");   // optimize girth
    H.display_stats();
	LDPC_Generator_Systematic G(&H);
    LDPC_Code C1(&H, &G);
    C1.save_code(FILENAME_IT);

  }

  if(0)
  {  // 16200 bits IRREGULAR (takes 3 hours to run)
	  cout << "========= IRREGULAR CODE 16200 BITS ==========" << endl;
	  LDPC_Parity_Irregular H;
	  H.generate(16200,
		  "0 0.21991 0.23328 0.02058 0 0.08543 0.06540 0.04767 0.01912"
		  "0 0 0 0 0 0 0 0 0 0.08064 0.22798",
		  "0 0 0 0 0 0 0 0.64854 0.34747 0.00399",
		  "rand",  // random unstructured matrix
		  "150 8"); // optimize
	  //LDPC_Code C(&H);
	  LDPC_Generator_Systematic G(&H);
	  LDPC_Code C(&H, &G);
	  C.save_code("../data/RU_16200.it");
  }

  ifstream  testfile;
  testfile.open( FILENAME_IT );
  if ( testfile == NULL )
  {
	  cout << "Done!" << endl << "Success to generate ldpc code file - \"random_3_6_16200.it\" in data path!" << endl ;
	  cout << "Please run dvb-s2.exe to demodulate and decode a data stream." << endl ;
	  return 0;
  }
  testfile.close();

  system( "pause" );

  return 0;

}
