/*!
 * \file
 * \brief Logarithmic and exponenential functions - header file
 * \author Tony Ottosson and Adam Piatyszek
 * 
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2005  (see AUTHORS file for a list of contributors)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 * -------------------------------------------------------------------------
 */

#ifndef LOGEXPFUNC_H
#define LOGEXPFUNC_H

#include <itpp/base/help_functions.h>

// Undefine log2 macro - IT++ has its own inline function 
#if defined (log2)
#undef log2
#endif


namespace itpp {

  //!\addtogroup logexpfunc
  //!@{

  // ----------------------------------------------------------------------
  // scalar functions
  // ----------------------------------------------------------------------

  //! Base-2 logarithm
  inline double log2(double x) 
  { 
    return (std::log(x) / 0.693147180559945309417); 
  }

  //! Base-b logarithm
  inline double logb(double b, double x) 
  { 
    return (std::log(x) / std::log(b)); 
  }

  //! Calculate two to the power of x (2^x); x is integer
  inline int pow2i(int x) { return ((x < 0) ? 0 : (1 << x)); }
  //! Calculate two to the power of x (2^x)
  inline double pow2(double x) { return pow(2.0, x); }

  //! Calculate ten to the power of x (10^x)
  inline double pow10(double x) { return pow(10.0, x); }

  //! Decibel of x (10*log10(x))
  inline double dB(double x) { return 10.0 * log10(x); }
  //! Inverse of decibel of x
  inline double inv_dB(double x) { return pow(10.0, 0.1 * x); }

  //! Calculate how many bits are needed to represent the integer n
  inline int needed_bits(int n)
  {
    int b = 0;
    it_assert(n > 0, "needed_bits(): n must be greater than zero!");
    n--; 
    while (n) { 
      n >>= 1; 
      b++; 
    }
    return b;
  }


  // ----------------------------------------------------------------------
  // functions on vectors and matrices
  // ----------------------------------------------------------------------

  //! Exp of the elements of a vector \c x
  inline vec exp(const vec &x) 
  { 
    return apply_function<double>(std::exp, x); 
  }
  //! Exp of the elements of a complex vector \c x
  inline cvec exp(const cvec &x) 
  {
    return apply_function<std::complex<double> >(std::exp, x);
  }
  //! Exp of the elements of a matrix \c m
  inline mat exp(const mat &m) 
  { 
    return apply_function<double>(std::exp, m); 
  }
  //! Exp of the elements of a complex matrix \c m
  inline cmat exp(const cmat &m) 
  { 
    return apply_function<std::complex<double> >(std::exp, m);
  }

  //! Calculates x to the power of y (x^y)
  inline vec pow(const double x, const vec &y)
  {
    return apply_function<double>(std::pow, x, y);
  }
  //! Calculates x to the power of y (x^y)
  inline mat pow(const double x, const mat &y)
  {
    return apply_function<double>(std::pow, x, y);
  }
  //! Calculates x to the power of y (x^y)
  inline vec pow(const vec &x, const double y)
  {
    return apply_function<double>(std::pow, x, y);
  }
  //! Calculates x to the power of y (x^y)
  inline mat pow(const mat &x, const double y)
  {
    return apply_function<double>(std::pow, x, y);
  }

  //! Calculates two to the power of x (2^x)
  inline vec pow2(const vec &x)
  {
    return apply_function<double>(pow2, x);
  }
  //! Calculates two to the power of x (2^x)
  inline mat pow2(const mat &x)
  {
    return apply_function<double>(pow2, x);
  }

  //! Calculates ten to the power of x (10^x)
  inline vec pow10(const vec &x)
  {
    return apply_function<double>(pow10, x);
  }
  //! Calculates ten to the power of x (10^x)
  inline mat pow10(const mat &x)
  {
    return apply_function<double>(pow10, x);
  }

  //! The natural logarithm of the elements
  inline vec log(const vec &x) 
  {
    return apply_function<double>(std::log, x);
  }
  //! The natural logarithm of the elements
  inline mat log(const mat &x) 
  {
    return apply_function<double>(std::log, x);
  }
  //! The natural logarithm of the elements
  inline cvec log(const cvec &x) 
  {
    return apply_function<std::complex<double> >(std::log, x);
  }
  //! The natural logarithm of the elements
  inline cmat log(const cmat &x) 
  {
    return apply_function<std::complex<double> >(std::log, x);
  }

  //! log-2 of the elements
  inline vec log2(const vec &x) 
  {
    return apply_function<double>(itpp::log2, x);
  }
  //! log-2 of the elements
  inline mat log2(const mat &x) 
  {
    return apply_function<double>(itpp::log2, x);
  }

  //! log-10 of the elements
  inline vec log10(const vec &x) 
  {
    return apply_function<double>(std::log10, x);
  }
  //! log-10 of the elements
  inline mat log10(const mat &x) 
  {
    return apply_function<double>(std::log10, x);
  }

  //! log-b of \c x
  inline vec logb(double b, const vec &x) 
  {
    return apply_function<double>(itpp::logb, b, x);
  }
  //! log-b of \c x
  inline mat logb(double b, const mat &x) 
  {
    return apply_function<double>(itpp::logb, b, x);
  }

  //! Calculates 10*log10(x)
  inline vec dB(const vec &x) 
  {
    return apply_function<double>(dB, x);
  }
  //! Calculates 10*log10(x)
  inline mat dB(const mat &x) 
  {
    return apply_function<double>(dB, x);
  }

  //! Calulates the inverse of dB, 10^(x/10)
  inline vec inv_dB(const vec &x) 
  {
    return apply_function<double>(inv_dB, x);
  }
  //! Calculates the inverse of dB, 10^(x/10)
  inline mat inv_dB(const mat &x) 
  {
    return apply_function<double>(inv_dB, x);
  }

  //! Calculate the numbers of bits needed for each symbol in a vector
  inline ivec needed_bits(const ivec& v) 
  {
    return apply_function<int>(needed_bits, v);
  }

  //!@}

} // namespace itpp

#endif // #ifndef LOGEXPFUNC_H



