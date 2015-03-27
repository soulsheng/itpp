/* -*- c++ -*- */
/* 
 * Copyright 2014 Ron Economos.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_DVBS2_BCH_BB_IMPL_H
#define INCLUDED_DVBS2_BCH_BB_IMPL_H

#include "bch_bb.h"

    class bch_bb_impl : public bch_bb
    {
     private:
      unsigned int kbch;
      unsigned int nbch;
      unsigned int bch_code;
      unsigned int m_poly_n_8[4];
      unsigned int m_poly_n_10[5];
      unsigned int m_poly_n_12[6];
      unsigned int m_poly_s_12[6];
      int poly_mult(const int*, int, const int*, int, int*);
      void poly_pack(const int*, unsigned int*, int);
      void poly_reverse(int*, int*, int);
      inline void reg_4_shift(unsigned int*);
      inline void reg_5_shift(unsigned int*);
      inline void reg_6_shift(unsigned int*);
      void bch_poly_build_tables(void);

     public:
      bch_bb_impl(CODE_RATE rate, FRAME_TYPE framesize);
      ~bch_bb_impl();

      // Where all the action really happens
      void forecast (int noutput_items, int* ninput_items_required);

      int general_work(int noutput_items,
		       int* ninput_items,
		       const void* input_items,
		       void* output_items);
    };

#endif /* INCLUDED_DVBS2_BCH_BB_IMPL_H */

