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


#ifndef INCLUDED_DVBS2_BBHEADER_BB_H
#define INCLUDED_DVBS2_BBHEADER_BB_H

#include "modulatorDefinition.h"

    /*!
     * \brief <+description of block+>
     * \ingroup dvbs2
     *
     */
    class bbheader_bb 
    {
     public:

      /*!
       * \brief Return a shared_ptr to a new instance of dvbs2::bbheader_bb.
       *
       * To avoid accidental use of raw pointers, dvbs2::bbheader_bb's
       * constructor is in a private implementation
       * class. dvbs2::bbheader_bb::make is the public interface for
       * creating new instances.
       */
      static bbheader_bb* make(CODE_RATE rate, Rolloff_Factor rolloff, FRAME_TYPE framesize);
	  virtual void general_work(int noutput_items,
		  const void *input_items,
		  void *output_items){ }

	  virtual bool verify( const void *input_items ){ return true;}

	  virtual int getkBCH()	{ return 0;}
	};

#endif /* INCLUDED_DVBS2_BBHEADER_BB_H */

