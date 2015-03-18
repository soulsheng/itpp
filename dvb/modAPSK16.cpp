/*!
 * \file
 * \brief One- and two-dimensional modulators - source file
 * \author Tony Ottosson and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#include "modAPSK16.h"
#include <itpp/comm/commfunc.h>
#include <itpp/base/math/elem_math.h>
#include <itpp/base/specmat.h>


namespace itpp
{

// ----------------------------------------------------------------------
// APSK16(4,12)
// ----------------------------------------------------------------------

void APSK16::set_M(int Mary)
{
  k = levels2bits(Mary);
  M = Mary;
  it_assert(pow2i(k) == M, "PSK::set_M(): M is not a power of 2");

  symbols.set_size(M);
  bitmap = graycode(k);

  float r1, r2;

  CODE_RATE rate = C3_4;

  FRAME_TYPE framesize = FECFRAME_NORMAL;
  r2 = 1.0f;
  if (framesize == FECFRAME_NORMAL)
  {
	  switch(rate)
	  {
	  case C2_3:
		  r1 = r2 / 3.15;
		  break;
	  case C3_4:
		  r1 = r2 / 2.85;
		  break;
	  case C4_5:
		  r1 = r2 / 2.75;
		  break;
	  case C5_6:
		  r1 = r2 / 2.70;
		  break;
	  case C8_9:
		  r1 = r2 / 2.60;
		  break;
	  case C9_10:
		  r1 = r2 / 2.57;
		  break;
	
	  default:
		  r1 = 0;
		  break;
	  }
  }
  else
  {
	  switch(rate)
	  {
	  case C2_3:
		  r1 = r2 / 3.15;
		  break;
	  case C3_4:
		  r1 = r2 / 2.85;
		  break;
	  case C4_5:
		  r1 = r2 / 2.75;
		  break;
	  case C5_6:
		  r1 = r2 / 2.70;
		  break;
	  case C8_9:
		  r1 = r2 / 2.60;
		  break;
	 
	  default:
		  r1 = 0;
		  break;
	  }
  }

  cvec& m_16apsk = symbols;
  double M_PI = pi;

  m_16apsk[0].real( r2 * cos(M_PI / 4.0));
  m_16apsk[0].imag( r2 * sin(M_PI / 4.0));
  m_16apsk[1].real( r2 * cos(-M_PI / 4.0));
  m_16apsk[1].imag( r2 * sin(-M_PI / 4.0));
  m_16apsk[2].real( r2 * cos(3 * M_PI / 4.0));
  m_16apsk[2].imag( r2 * sin(3 * M_PI / 4.0));
  m_16apsk[3].real( r2 * cos(-3 * M_PI / 4.0));
  m_16apsk[3].imag( r2 * sin(-3 * M_PI / 4.0));
  m_16apsk[4].real( r2 * cos(M_PI / 12.0));
  m_16apsk[4].imag( r2 * sin(M_PI / 12.0));
  m_16apsk[5].real( r2 * cos(-M_PI / 12.0));
  m_16apsk[5].imag( r2 * sin(-M_PI / 12.0));
  m_16apsk[6].real( r2 * cos(11 * M_PI / 12.0));
  m_16apsk[6].imag( r2 * sin(11 * M_PI / 12.0));
  m_16apsk[7].real( r2 * cos(-11 * M_PI / 12.0));
  m_16apsk[7].imag( r2 * sin(-11 * M_PI / 12.0));
  m_16apsk[8].real( r2 * cos(5 * M_PI / 12.0));
  m_16apsk[8].imag( r2 * sin(5 * M_PI / 12.0));
  m_16apsk[9].real( r2 * cos(-5 * M_PI / 12.0));
  m_16apsk[9].imag( r2 * sin(-5 * M_PI / 12.0));
  m_16apsk[10].real( r2 * cos(7 * M_PI / 12.0));
  m_16apsk[10].imag( r2 * sin(7 * M_PI / 12.0));
  m_16apsk[11].real( r2 * cos(-7 * M_PI / 12.0));
  m_16apsk[11].imag( r2 * sin(-7 * M_PI / 12.0));
  m_16apsk[12].real( r1 * cos(M_PI / 4.0));
  m_16apsk[12].imag( r1 * sin(M_PI / 4.0));
  m_16apsk[13].real( r1 * cos(-M_PI / 4.0));
  m_16apsk[13].imag( r1 * sin(-M_PI / 4.0));
  m_16apsk[14].real( r1 * cos(3 * M_PI / 4.0));
  m_16apsk[14].imag( r1 * sin(3 * M_PI / 4.0));
  m_16apsk[15].real( r1 * cos(-3 * M_PI / 4.0));
  m_16apsk[15].imag( r1 * sin(-3 * M_PI / 4.0));

  // W1, W2
  double rightA = r2 * cos( M_PI/4 );
  double leftA2 = r2 * cos( M_PI*5/12 );
  double leftA1 = r1 * cos( M_PI/4 );
  double leftA = leftA2 > leftA1 ? leftA2 : leftA1;
  W = ( rightA + leftA )*0.5f ;

  // factor
  double alpha = 1.2703f;
  double amp = (4*r1 + 12*r2)/16;
  factor = alpha * amp ;


  calculate_softbit_matrices();

  setup_done = true;
}

void APSK16::demodulate_bits(const cvec &signal, bvec &out) const
{
  it_assert_debug(setup_done, "PSK::demodulate_bits(): Modulator not ready.");
  int est_symbol;
  double ang, temp;

  out.set_size(k*signal.size(), false);

  for (int i = 0; i < signal.size(); i++) {
   
    est_symbol = 0;

	double distance2 = 1000.0f;
	for(int j=0;j<M;j++)
	{
		double distTemp = std::norm( signal(i) - symbols[j] );
		if( distTemp<distance2 )
		{
			distance2 = distTemp;
			est_symbol = j;
		}
	}

    out.replace_mid(i*k, dec2bin(k,est_symbol));
  }
}

bvec APSK16::demodulate_bits(const cvec &signal) const
{
  bvec out;
  demodulate_bits(signal, out);
  return out;
}

void APSK16::demodulate_soft_bits(const cvec &rx_symbols, double N0,
                                vec &soft_bits, Soft_Method) const
{
  soft_bits.set_size(k * rx_symbols.size());
 
  for (int i = 0; i < rx_symbols.size(); i++) {
	double a = std::real( rx_symbols(i) );
	double b = std::imag( rx_symbols(i) );
   
	soft_bits((i * k) + 3) = W - fabs(a);
	soft_bits((i * k) + 2) = W - fabs(b);
    soft_bits((i * k) + 1) = -a;
    soft_bits(i * k) = -b;
  }
}

vec APSK16::demodulate_soft_bits(const cvec &rx_symbols, double N0,
                               Soft_Method method) const
{
  vec out;
  demodulate_soft_bits(rx_symbols, N0, out, method);
  return out;
}

void APSK16::modulate_bits( const bvec& bits, cvec& output ) const
{
	int no_symbols = bits.length() / k;
	output.set_size(no_symbols);
	for (int i = 0; i < no_symbols; i++) {
		bvec bits_k = bits.mid(i * k, k);
		unsigned int symbol_in = bin2dec(bits_k);
		output(i) = symbols( symbol_in );
	}
}

cvec APSK16::modulate_bits( const bvec& bits ) const
{
	cvec output;
	modulate_bits(bits, output);
	return output;
}


} // namespace itpp
