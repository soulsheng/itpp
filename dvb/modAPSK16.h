/*!
 * \file
 * \brief One- and two-dimensional modulators - header file
 * \author Tony Ottosson and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2011  (see AUTHORS file for a list of contributors)
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

#ifndef MODAPSK16_H
#define MODAPSK16_H

#include <itpp/base/mat.h>
#include <itpp/base/math/elem_math.h>
#include <itpp/base/math/log_exp.h>
#include <itpp/base/converters.h>
#include <itpp/base/math/min_max.h>

#include <itpp/comm/modulator.h>

namespace itpp
{

// ----------------------------------------------------------------------
// PSK : Modulator<std::complex<double> >
// ----------------------------------------------------------------------

/*!
  \ingroup modulators
  \brief M-ary APSK modulator.

  This class implements the M-ary PSK modulator with \f$M = 2^k\f$
  constellation points, where \f$k = 1, 2, \ldots \f$. The symbol
  numbering is counter clockwise starting from the real axis, i.e. symbol
  \f$(1, 0)\f$. The bitmap is Gray encoded. The symbol energy is
  normalized to 1.

  Beside hard demapping, this class can also perform soft demodulation,
  calculating the log-MAP estimate of the individual bits. To use it
  properly the received symbols should be equal to: \f[r_k = c_k s_k +
  n_k,\f] where \f$c_k\f$ is the real or complex channel gain, \f$s_k\f$
  is the transmitted constellation symbol, and \f$n_k\f$ is the AWGN of
  the channel (with variance \f$N_0\f$).

  It is also assumed that the channel estimates are perfect when
  calculating the soft bits.
*/
class APSK16 : public Modulator<std::complex<double> >
{
public:
  //! Default Constructor
  APSK16() {}
  //! Class constructor
  APSK16(int M) { set_M(M); }
  //! Destructor
  virtual ~APSK16() { }
  //! Change the size of the signal constellation
  void set_M(int M);


  //! Modulation of bits
  virtual void modulate_bits(const bvec& bits, cvec& output) const;
  //! Modulation of bits
  virtual cvec modulate_bits(const bvec& bits) const;

  //! Hard demodulation of bits
  void demodulate_bits(const cvec& signal, bvec& bits) const;
  //! Hard demodulation of bits
  bvec demodulate_bits(const cvec& signal) const;

  /*!
    \brief Soft demodulator for AWGN channel

    This function calculates the log-MAP estimate assuming equally likely
    bits transmitted: \f[\log \left( \frac{P(b=0|r)}{P(b=1|r)} \right) =
    \frac{2 \sqrt{2}}{N_0} \Im\{r_k \exp \left(j \frac{\Pi}{4} \right)
    \}\f] and \f[\log \left( \frac{P(b=0|r)}{P(b=1|r)} \right) = \frac{2
    \sqrt{2}}{N_0} \Re\{r_k \exp \left(j \frac{\Pi}{4} \right) \}\f]
    depending on the bit positon in the QPSK symbol.

    \param rx_symbols The received noisy constellation symbols, \f$r\f$
    \param N0 The spectral density of the AWGN noise, \f$n\f$
    \param soft_bits The soft bits calculated using the expression above
    \param method The method used for demodulation (LOGMAP or APPROX)

    \note For soft demodulation it is suggested to use the
    N-dimensional modulator (Modulator_ND) instead, which is based
    on the QLLR (quantized) arithmetic and therefore is
    faster. Please note, however, that mixed use of \c
    Modulator_1D/\c Modulator_2D and \c Modulator_ND is not advised.
  */
  virtual void demodulate_soft_bits(const cvec& rx_symbols, double N0,
                                    vec& soft_bits,
                                    Soft_Method method = LOGMAP) const;
  //! Soft demodulator for AWGN channel
  vec demodulate_soft_bits(const cvec& rx_symbols, double N0,
                           Soft_Method method = LOGMAP) const;

  double	W ; // W1 = W2
  double	factor;
};

} // namespace itpp

#endif // #ifndef MODULATOR_H
