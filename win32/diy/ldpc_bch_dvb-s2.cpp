#include <itpp/itcomm.h>
#include <sstream>

using namespace std;
using namespace itpp;

#define		FILENAME_IT		"RU_16200.it"
#define		EBNO			1.8
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

bool syndrome_check(const QLLRvec &LLR,
	int ncheck, 
	ivec& sumX2, 
	ivec& V) ;

int bp_decode(const QLLRvec &LLRin, int *LLRout,
	int nvar, int ncheck, 
	ivec& C, ivec& V, ivec& sumX1, ivec& sumX2, ivec& iind, ivec& jind,	// Parity check matrix parameterization
	QLLRvec& mvc, QLLRvec&mcv,	// temporary storage for decoder (memory allocated when codec defined)
	//LLR_calc_unit& llrcalc,		//!< LLR calculation unit
	short int Dint1, short int Dint2, short int Dint3,	//! Decoder (lookup-table) parameters
	ivec& logexp_table,		//! The lookup tables for the decoder
	bool psc = true,			//!< check syndrom after each iteration
	int max_iters = 50 );		//!< Maximum number of iterations

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

	Real_Timer	timer, timerStep;
	vec			timerValue(COUNT_REPEAT);

	vec			timerStepValue(COUNT_REPEAT);
	ivec		countIteration(COUNT_REPEAT);

	int nldpc = ldpc.get_nvar();
	int kldpc = ldpc.get_ninfo();             // number of bits per codeword

	int nSplit = kldpc / N_BCH;
	int Nbch = nSplit * N_BCH;

	int Kbch = nSplit * K_BCH;

	ferc.set_blocksize(Kbch);

	int*	pIntLLROut = (int*)malloc( nldpc*sizeof(int) );

    for (int64_t i = 0; i < COUNT_REPEAT; i ++) 
	{
		timer.reset();
		timer.start();

		// step 1: input message
		bvec bitsinBCHEnc = randb(Kbch);

		// step 2: bch encode
		bvec bitsoutBCHEnc = zeros_b(Nbch);
		for (int j = 0; j < nSplit; j++)
		{
			bvec bitsinBCHOne = bitsinBCHEnc.mid(j*K_BCH, K_BCH);
			bvec encodedBCHOne = bch.encode(bitsinBCHOne);
			bitsoutBCHEnc.set_subvector(j*N_BCH, encodedBCHOne);
		}

		bvec bitsinLDPCEnc = zeros_b(kldpc);
		bitsinLDPCEnc.set_subvector(0, bitsoutBCHEnc);


		// step 3: ldpc encode
		bvec bitsoutLDPCEnc = ldpc.encode(bitsinLDPCEnc);
		
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




		// step 7: ldpc Decode the received bits
		QLLRvec llr(nldpc);
		QLLRvec llrIn = ldpc.get_llrcalc().to_qllr(softbits);
		memset( pIntLLROut, 0, sizeof(int)*nldpc );

		timerStep.reset();
		timerStep.start();

		countIteration[i] = bp_decode(llrIn, pIntLLROut,
			ldpc.nvar, ldpc.ncheck, 
			ldpc.C, ldpc.V, ldpc.sumX1, ldpc.sumX2, ldpc.iind, ldpc.jind,	// Parity check matrix parameterization
			ldpc.mvc, ldpc.mcv,	// temporary storage for decoder (memory allocated when codec defined)
			//ldpc.llrcalc );		//!< LLR calculation unit
			ldpc.llrcalc.Dint1, ldpc.llrcalc.Dint2, ldpc.llrcalc.Dint3,	//! Decoder (lookup-table) parameters
			ldpc.llrcalc.logexp_table);		//! The lookup tables for the decoder

		timerStep.stop();
		timerStepValue[i] = timerStep.get_time() ;

		for (int i=0;i<nldpc;i++)
			llr[i] = pIntLLROut[i];

		bvec bitsoutLDPCDec = llr < 0;
		bvec bitsinBCHDec = bitsoutLDPCDec.left(Nbch);
		//      bvec bitsout = C.decode(softbits); // (only systematic bits)


		// step 8: bch decode
		bvec bitsoutBCHDec = zeros_b(Kbch);
		for (int j = 0; j < nSplit; j++)
		{
			bvec bitsinBCHDecOne = bitsinBCHDec.mid(j*N_BCH, N_BCH);
			bvec bitsoutBCHDecOne = bch.decode(bitsinBCHDecOne);
			bitsoutBCHDec.set_subvector(j*K_BCH, bitsoutBCHDecOne);
		}

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

	double timerAverageAll = 0.0f, timerStepAverage = 0.0f;
	for (int i=0;i<COUNT_REPEAT;i++)
	{
		cout << timerValue[i] << " s, " ;
		timerAverageAll += timerValue[i];
	}
	cout << endl << endl ;

	for (int i=0;i<COUNT_REPEAT;i++)
	{
		cout << timerStepValue[i] << " s, " ;
		timerStepAverage += timerStepValue[i];
	}
	cout << endl << endl ;

	double countIterationAverage = 0.0f;
	for (int i=0;i<COUNT_REPEAT;i++)
	{
		cout << countIteration[i] << " iteration, " ;
		countIterationAverage += countIteration[i];
	}
	cout << endl << endl ;

	cout << endl << timerAverageAll/COUNT_REPEAT << " s Average all" << endl ;
	
	cout << endl << timerStepAverage/COUNT_REPEAT << " s Average step decode ldpc" << endl ;
	
	cout << endl << countIterationAverage/COUNT_REPEAT << " iteration Average in decode ldpc" << endl ;

	free( pIntLLROut );

	return 0;
}

bool syndrome_check(int *LLR,
	int ncheck, 
	ivec& sumX2, 
	ivec& V ) 
{
	// Please note the IT++ convention that a sure zero corresponds to
	// LLR=+infinity
	int i, j, synd, vi;

	for (j = 0; j < ncheck; j++) {
		synd = 0;
		int vind = j; // tracks j+i*ncheck
		for (i = 0; i < sumX2(j); i++) {
			vi = V(vind);
			if (LLR[vi] < 0) {
				synd++;
			}
			vind += ncheck;
		}
		if ((synd&1) == 1) {
			return false;  // codeword is invalid
		}
	}
	return true;   // codeword is valid
}

QLLR  logexp(QLLR x,
	short int Dint1, short int Dint2, short int Dint3,	//! Decoder (lookup-table) parameters
	ivec& logexp_table )		//! The lookup tables for the decoder
{
	int ind = x >> Dint3;
	if (ind >= Dint2) // outside table
		return 0;

	// Without interpolation
	return logexp_table(ind);
}

QLLR Boxplus(QLLR a, QLLR b,
	short int Dint1, short int Dint2, short int Dint3,	//! Decoder (lookup-table) parameters
	ivec& logexp_table )		//! The lookup tables for the decoder
{
	QLLR a_abs = (a > 0 ? a : -a);
	QLLR b_abs = (b > 0 ? b : -b);
	QLLR minabs = (a_abs > b_abs ? b_abs : a_abs);
	QLLR term1 = (a > 0 ? (b > 0 ? minabs : -minabs)
		: (b > 0 ? -minabs : minabs));

	if (Dint2 == 0) {  // logmax approximation - avoid looking into empty table
		// Don't abort when overflowing, just saturate the QLLR
		if (term1 > QLLR_MAX) {
			return QLLR_MAX;
		}
		if (term1 < -QLLR_MAX) {
			return -QLLR_MAX;
		}
		return term1;
	}

	QLLR apb = a + b;
	QLLR term2 = logexp((apb > 0 ? apb : -apb), Dint1, Dint2, Dint3, logexp_table);
	QLLR amb = a - b;
	QLLR term3 = logexp((amb > 0 ? amb : -amb), Dint1, Dint2, Dint3, logexp_table);
	QLLR result = term1 + term2 - term3;

	// Don't abort when overflowing, just saturate the QLLR
	if (result > QLLR_MAX) {
		return QLLR_MAX;
	}
	if (result < -QLLR_MAX) {
		return -QLLR_MAX;
	}
	return result;
}

int bp_decode(const QLLRvec &LLRin, int *LLRout,
	int nvar, int ncheck, 
	ivec& C, ivec& V, ivec& sumX1, ivec& sumX2, ivec& iind, ivec& jind,	// Parity check matrix parameterization
	QLLRvec& mvc, QLLRvec&mcv,	// temporary storage for decoder (memory allocated when codec defined)
	//LLR_calc_unit& llrcalc,		//!< LLR calculation unit
	short int Dint1, short int Dint2, short int Dint3,	//! Decoder (lookup-table) parameters
	ivec& logexp_table,		//! The lookup tables for the decoder
	bool psc /*= true*/,			//!< check syndrom after each iteration
	int max_iters /*= 50*/ )		//!< Maximum number of iterations
{
  // Note the IT++ convention that a sure zero corresponds to
  // LLR=+infinity

  //LLRout.set_size(LLRin.size());

  //! Maximum check node degree that the class can handle
  static const int max_cnd = 200;

  // allocate temporary variables used for the check node update
  ivec jj(max_cnd);
  QLLRvec m(max_cnd);
  QLLRvec ml(max_cnd);
  QLLRvec mr(max_cnd);
  
  // initial step
  for (int i = 0; i < nvar; i++) {
    int index = i;
    for (int j = 0; j < sumX1(i); j++) {
      mvc[index] = LLRin(i);
      index += nvar;
    }
  }

  bool is_valid_codeword = false;
  int iter = 0;
  do {
    iter++;
    if (nvar >= 100000) { it_info_no_endl_debug("."); }
    // --------- Step 1: check to variable nodes ----------
    for (int j = 0; j < ncheck; j++) {
      // The check node update calculations are hardcoded for degrees
      // up to 6.  For larger degrees, a general algorithm is used.
      switch (sumX2(j)) {
      case 0:
        it_error("LDPC_Code::bp_decode(): sumX2(j)=0");
      case 1:
        it_error("LDPC_Code::bp_decode(): sumX2(j)=1");
      case 2: {
        mcv[j+ncheck] = mvc[jind[j]];
        mcv[j] = mvc[jind[j+ncheck]];
        break;
      }
      case 3: {
        int j0 = j;
        QLLR m0 = mvc[jind[j0]];
        int j1 = j0 + ncheck;
        QLLR m1 = mvc[jind[j1]];
        int j2 = j1 + ncheck;
        QLLR m2 = mvc[jind[j2]];
        mcv[j0] = Boxplus(m1, m2, Dint1, Dint2, Dint3, logexp_table);
        mcv[j1] = Boxplus(m0, m2, Dint1, Dint2, Dint3, logexp_table);
        mcv[j2] = Boxplus(m0, m1, Dint1, Dint2, Dint3, logexp_table);
        break;
      }
      case 4: {
        int j0 = j;
        QLLR m0 = mvc[jind[j0]];
        int j1 = j0 + ncheck;
        QLLR m1 = mvc[jind[j1]];
        int j2 = j1 + ncheck;
        QLLR m2 = mvc[jind[j2]];
        int j3 = j2 + ncheck;
        QLLR m3 = mvc[jind[j3]];
        QLLR m01 = Boxplus(m0, m1, Dint1, Dint2, Dint3, logexp_table);
        QLLR m23 = Boxplus(m2, m3, Dint1, Dint2, Dint3, logexp_table);
        mcv[j0] = Boxplus(m1, m23, Dint1, Dint2, Dint3, logexp_table);
        mcv[j1] = Boxplus(m0, m23, Dint1, Dint2, Dint3, logexp_table);
        mcv[j2] = Boxplus(m01, m3, Dint1, Dint2, Dint3, logexp_table);
        mcv[j3] = Boxplus(m01, m2, Dint1, Dint2, Dint3, logexp_table);
        break;
      }
      case 5: {
        int j0 = j;
        QLLR m0 = mvc[jind[j0]];
        int j1 = j0 + ncheck;
        QLLR m1 = mvc[jind[j1]];
        int j2 = j1 + ncheck;
        QLLR m2 = mvc[jind[j2]];
        int j3 = j2 + ncheck;
        QLLR m3 = mvc[jind[j3]];
        int j4 = j3 + ncheck;
        QLLR m4 = mvc[jind[j4]];
        QLLR m01 = Boxplus(m0, m1, Dint1, Dint2, Dint3, logexp_table);
        QLLR m02 = Boxplus(m01, m2, Dint1, Dint2, Dint3, logexp_table);
        QLLR m34 = Boxplus(m3, m4, Dint1, Dint2, Dint3, logexp_table);
        QLLR m24 = Boxplus(m2, m34, Dint1, Dint2, Dint3, logexp_table);
        mcv[j0] = Boxplus(m1, m24, Dint1, Dint2, Dint3, logexp_table);
        mcv[j1] = Boxplus(m0, m24, Dint1, Dint2, Dint3, logexp_table);
        mcv[j2] = Boxplus(m01, m34, Dint1, Dint2, Dint3, logexp_table);
        mcv[j3] = Boxplus(m02, m4, Dint1, Dint2, Dint3, logexp_table);
        mcv[j4] = Boxplus(m02, m3, Dint1, Dint2, Dint3, logexp_table);
        break;
      }
      case 6: {
        int j0 = j;
        QLLR m0 = mvc[jind[j0]];
        int j1 = j0 + ncheck;
        QLLR m1 = mvc[jind[j1]];
        int j2 = j1 + ncheck;
        QLLR m2 = mvc[jind[j2]];
        int j3 = j2 + ncheck;
        QLLR m3 = mvc[jind[j3]];
        int j4 = j3 + ncheck;
        QLLR m4 = mvc[jind[j4]];
        int j5 = j4 + ncheck;
        QLLR m5 = mvc[jind[j5]];
        QLLR m01 = Boxplus(m0, m1, Dint1, Dint2, Dint3, logexp_table);
        QLLR m23 = Boxplus(m2, m3, Dint1, Dint2, Dint3, logexp_table);
        QLLR m45 = Boxplus(m4, m5, Dint1, Dint2, Dint3, logexp_table);
        QLLR m03 = Boxplus(m01, m23, Dint1, Dint2, Dint3, logexp_table);
        QLLR m25 = Boxplus(m23, m45, Dint1, Dint2, Dint3, logexp_table);
        QLLR m0145 = Boxplus(m01, m45, Dint1, Dint2, Dint3, logexp_table);
        mcv[j0] = Boxplus(m1, m25, Dint1, Dint2, Dint3, logexp_table);
        mcv[j1] = Boxplus(m0, m25, Dint1, Dint2, Dint3, logexp_table);
        mcv[j2] = Boxplus(m0145, m3, Dint1, Dint2, Dint3, logexp_table);
        mcv[j3] = Boxplus(m0145, m2, Dint1, Dint2, Dint3, logexp_table);
        mcv[j4] = Boxplus(m03, m5, Dint1, Dint2, Dint3, logexp_table);
        mcv[j5] = Boxplus(m03, m4, Dint1, Dint2, Dint3, logexp_table);
        break;
      }
      default: {
        int nodes = sumX2(j);
        if( nodes > max_cnd ) {
          std::ostringstream m_sout;
          m_sout << "check node degrees >" << max_cnd << " not supported in this version";
          it_error( m_sout.str() );
        }

        nodes--;
        jj[0] = j;
        m[0] = mvc[jind[jj[0]]];
        for(int i = 1; i <= nodes; i++ ) {
          jj[i] = jj[i-1] + ncheck;
          m[i] = mvc[jind[jj[i]]];
        }

	// compute partial sums from the left and from the right
        ml[0] = m[0];
        mr[0] = m[nodes];
        for(int i = 1; i < nodes; i++ ) {
          ml[i] = Boxplus( ml[i-1], m[i], Dint1, Dint2, Dint3, logexp_table );
          mr[i] = Boxplus( mr[i-1], m[nodes-i], Dint1, Dint2, Dint3, logexp_table );
        }

	// merge partial sums
        mcv[jj[0]] = mr[nodes-1];
        mcv[jj[nodes]] = ml[nodes-1];
        for(int i = 1; i < nodes; i++ )
          mcv[jj[i]] = Boxplus( ml[i-1], mr[nodes-1-i], Dint1, Dint2, Dint3, logexp_table );
      }
      }  // switch statement
    }
    
    // step 2: variable to check nodes
    for (int i = 0; i < nvar; i++) {
      switch (sumX1(i)) {
      case 0:
        it_error("LDPC_Code::bp_decode(): sumX1(i)=0");
      case 1: {
        /* This case is rare but apparently occurs for codes used in
	   the DVB-T2 standard.
	 */
	QLLR m0 = mcv[iind[i]];
        mvc[i] = LLRin(i);
        LLRout[i] = LLRin(i) + m0;
        break;
      }
      case 2: {
        QLLR m0 = mcv[iind[i]];
        int i1 = i + nvar;
        QLLR m1 = mcv[iind[i1]];
        mvc[i] = LLRin(i) + m1;
        mvc[i1] = LLRin(i) + m0;
        LLRout[i] = mvc[i1] + m1;
        break;
      }
      case 3: {
        int i0 = i;
        QLLR m0 = mcv[iind[i0]];
        int i1 = i0 + nvar;
        QLLR m1 = mcv[iind[i1]];
        int i2 = i1 + nvar;
        QLLR m2 = mcv[iind[i2]];
        LLRout[i] = LLRin(i) + m0 + m1 + m2;
        mvc[i0] = LLRout[i] - m0;
        mvc[i1] = LLRout[i] - m1;
        mvc[i2] = LLRout[i] - m2;
        break;
      }
      case 4: {
        int i0 = i;
        QLLR m0 = mcv[iind[i0]];
        int i1 = i0 + nvar;
        QLLR m1 = mcv[iind[i1]];
        int i2 = i1 + nvar;
        QLLR m2 = mcv[iind[i2]];
        int i3 = i2 + nvar;
        QLLR m3 = mcv[iind[i3]];
        LLRout[i] = LLRin(i) + m0 + m1 + m2 + m3;
        mvc[i0] = LLRout[i] - m0;
        mvc[i1] = LLRout[i] - m1;
        mvc[i2] = LLRout[i] - m2;
        mvc[i3] = LLRout[i] - m3;
        break;
      }
      default:   { // differential update
        QLLR mvc_temp = LLRin(i);
        int index_iind = i; // tracks i+jp*nvar
        for (int jp = 0; jp < sumX1(i); jp++) {
          mvc_temp +=  mcv[iind[index_iind]];
          index_iind += nvar;
        }
        LLRout[i] = mvc_temp;
        index_iind = i;  // tracks i+j*nvar
        for (int j = 0; j < sumX1[i]; j++) {
          mvc[index_iind] = mvc_temp - mcv[iind[index_iind]];
          index_iind += nvar;
        }
      }
      }
    }

	if (psc && syndrome_check(LLRout, ncheck, sumX2, V)) {
      is_valid_codeword = true;
      break;
    }
  }
  while (iter < max_iters);

  return (is_valid_codeword ? iter : -iter);
}
