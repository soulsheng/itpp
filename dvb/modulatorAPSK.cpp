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

#include "modulatorAPSK.h"
#include <itpp/comm/commfunc.h>
#include <itpp/base/math/elem_math.h>
#include <itpp/base/specmat.h>


namespace itpp
{

// ----------------------------------------------------------------------
// APSK32(4,12,16)
// ----------------------------------------------------------------------

enum CODE_RATE
{
	C3_4,
	C4_5,
	C5_6,
	C8_9,
	C9_10
};

void APSK32::set_M(int Mary)
{
  k = levels2bits(Mary);
  M = Mary;
  it_assert(pow2i(k) == M, "PSK::set_M(): M is not a power of 2");

  symbols.set_size(M);
  bitmap = graycode(k);

  float r1, r2;
  float r3 = 1.0f;

  CODE_RATE rate = C3_4;

  switch(rate)
  {
  case C3_4:
	  r1 = r3 / 5.27;
	  r2 = r1 * 2.84;
	  break;
  case C4_5:
	  r1 = r3 / 4.87;
	  r2 = r1 * 2.72;
	  break;
  case C5_6:
	  r1 = r3 / 4.64;
	  r2 = r1 * 2.64;
	  break;
  case C8_9:
	  r1 = r3 / 4.33;
	  r2 = r1 * 2.54;
	  break;
  case C9_10:
	  r1 = r3 / 4.30;
	  r2 = r1 * 2.53;
	  break;
  default:
	  r1 = 0;
	  r2 = 0;
	  break;
  }

  cvec& m_32apsk = symbols;
  double M_PI = pi;
  symbols[0].real( r2 * cos(M_PI / 4.0) );
  symbols[0].imag( r2 * sin(M_PI / 4.0));
  m_32apsk[1].real( r2 * cos(5 * M_PI / 12.0));
  m_32apsk[1].imag( r2 * sin(5 * M_PI / 12.0));
  m_32apsk[2].real( r2 * cos(-M_PI / 4.0));
  m_32apsk[2].imag( r2 * sin(-M_PI / 4.0));
  m_32apsk[3].real( r2 * cos(-5 * M_PI / 12.0));
  m_32apsk[3].imag( r2 * sin(-5 * M_PI / 12.0));
  m_32apsk[4].real( r2 * cos(3 * M_PI / 4.0));
  m_32apsk[4].imag( r2 * sin(3 * M_PI / 4.0));
  m_32apsk[5].real( r2 * cos(7 * M_PI / 12.0));
  m_32apsk[5].imag( r2 * sin(7 * M_PI / 12.0));
  m_32apsk[6].real( r2 * cos(-3 * M_PI / 4.0));
  m_32apsk[6].imag( r2 * sin(-3 * M_PI / 4.0));
  m_32apsk[7].real( r2 * cos(-7 * M_PI / 12.0));
  m_32apsk[7].imag( r2 * sin(-7 * M_PI / 12.0));
  m_32apsk[8].real( r3 * cos(M_PI / 8.0));
  m_32apsk[8].imag( r3 * sin(M_PI / 8.0));
  m_32apsk[9].real( r3 * cos(3 * M_PI / 8.0));
  m_32apsk[9].imag( r3 * sin(3 * M_PI / 8.0));
  m_32apsk[10].real( r3 * cos(-M_PI / 4.0));
  m_32apsk[10].imag( r3 * sin(-M_PI / 4.0));
  m_32apsk[11].real( r3 * cos(-M_PI / 2.0));
  m_32apsk[11].imag( r3 * sin(-M_PI / 2.0));
  m_32apsk[12].real( r3 * cos(3 * M_PI / 4.0));
  m_32apsk[12].imag( r3 * sin(3 * M_PI / 4.0));
  m_32apsk[13].real( r3 * cos(M_PI / 2.0));
  m_32apsk[13].imag( r3 * sin(M_PI / 2.0));
  m_32apsk[14].real( r3 * cos(-7 * M_PI / 8.0));
  m_32apsk[14].imag( r3 * sin(-7 * M_PI / 8.0));
  m_32apsk[15].real( r3 * cos(-5 * M_PI / 8.0));
  m_32apsk[15].imag( r3 * sin(-5 * M_PI / 8.0));
  m_32apsk[16].real( r2 * cos(M_PI / 12.0));
  m_32apsk[16].imag( r2 * sin(M_PI / 12.0));
  m_32apsk[17].real( r1 * cos(M_PI / 4.0));
  m_32apsk[17].imag( r1 * sin(M_PI / 4.0));
  m_32apsk[18].real( r2 * cos(-M_PI / 12.0));
  m_32apsk[18].imag( r2 * sin(-M_PI / 12.0));
  m_32apsk[19].real( r1 * cos(-M_PI / 4.0));
  m_32apsk[19].imag( r1 * sin(-M_PI / 4.0));
  m_32apsk[20].real( r2 * cos(11 * M_PI / 12.0));
  m_32apsk[20].imag( r2 * sin(11 * M_PI / 12.0));
  m_32apsk[21].real( r1 * cos(3 * M_PI / 4.0));
  m_32apsk[21].imag( r1 * sin(3 * M_PI / 4.0));
  m_32apsk[22].real( r2 * cos(-11 * M_PI / 12.0));
  m_32apsk[22].imag( r2 * sin(-11 * M_PI / 12.0));
  m_32apsk[23].real( r1 * cos(-3 * M_PI / 4.0));
  m_32apsk[23].imag( r1 * sin(-3 * M_PI / 4.0));
  m_32apsk[24].real( r3 * cos(0.0));
  m_32apsk[24].imag( r3 * sin(0.0));
  m_32apsk[25].real( r3 * cos(M_PI / 4.0));
  m_32apsk[25].imag( r3 * sin(M_PI / 4.0));
  m_32apsk[26].real( r3 * cos(-M_PI / 8.0));
  m_32apsk[26].imag( r3 * sin(-M_PI / 8.0));
  m_32apsk[27].real( r3 * cos(-3 * M_PI / 8.0));
  m_32apsk[27].imag( r3 * sin(-3 * M_PI / 8.0));
  m_32apsk[28].real( r3 * cos(7 * M_PI / 8.0));
  m_32apsk[28].imag( r3 * sin(7 * M_PI / 8.0));
  m_32apsk[29].real( r3 * cos(5 * M_PI / 8.0));
  m_32apsk[29].imag( r3 * sin(5 * M_PI / 8.0));
  m_32apsk[30].real( r3 * cos(M_PI));
  m_32apsk[30].imag( r3 * sin(M_PI));
  m_32apsk[31].real( r3 * cos(-3 * M_PI / 4.0));
  m_32apsk[31].imag( r3 * sin(-3 * M_PI / 4.0));

  calculate_softbit_matrices();

  setup_done = true;
}


void APSK32::demodulate_bits(const cvec &signal, bvec &out) const
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

    out.replace_mid(i*k, dec2bin(k, est_symbol));
  }
}

bvec APSK32::demodulate_bits(const cvec &signal) const
{
  bvec out;
  demodulate_bits(signal, out);
  return out;
}


// ----------------------------------------------------------------------
// APSK32
// ----------------------------------------------------------------------

void APSK32::demodulate_soft_bits(const cvec &rx_symbols, double N0,
                                vec &soft_bits, Soft_Method) const
{
  soft_bits.set_size(k * rx_symbols.size());
  std::complex<double> temp;
  double factor = 2 * std::sqrt(2.0) / N0;
  std::complex<double> exp_pi4 = std::complex<double>(std::cos(pi / 4),
                                 std::sin(pi / 4));
  for (int i = 0; i < rx_symbols.size(); i++) {
    temp = rx_symbols(i) * exp_pi4;
    soft_bits((i << 1) + 1) = std::real(temp) * factor;
    soft_bits(i << 1) = std::imag(temp) * factor;
  }
}

vec APSK32::demodulate_soft_bits(const cvec &rx_symbols, double N0,
                               Soft_Method method) const
{
  vec out;
  demodulate_soft_bits(rx_symbols, N0, out, method);
  return out;
}


void APSK32::modulate_bits( const bvec& bits, cvec& output ) const
{
	int no_symbols = bits.length() / k;
	output.set_size(no_symbols);
	for (int i = 0; i < no_symbols; i++) {
		bvec bits_k = bits.mid(i * k, k);
		unsigned int symbol_in = bin2dec(bits_k);
		output(i) = symbols( symbol_in );
	}
}

cvec APSK32::modulate_bits( const bvec& bits ) const
{
	cvec output;
	modulate_bits(bits, output);
	return output;
}


} // namespace itpp
