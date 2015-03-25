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

#include "bbscrambler_bb_impl.h"
#include <stdio.h>


    bbscrambler_bb*
    bbscrambler_bb::make(CODE_RATE rate, FRAME_TYPE framesize)
    {
      return new bbscrambler_bb_impl(rate, framesize);
    }

    /*
     * The private constructor
     */
    bbscrambler_bb_impl::bbscrambler_bb_impl(CODE_RATE rate, FRAME_TYPE framesize)
    {
        if (framesize == FECFRAME_NORMAL)
        {
            switch (rate)
            {
                case C1_4:
                    kbch = 16008;
                    break;
                case C1_3:
                    kbch = 21408;
                    break;
                case C2_5:
                    kbch = 25728;
                    break;
                case C1_2:
                    kbch = 32208;
                    break;
                case C3_5:
                    kbch = 38688;
                    break;
                case C2_3:
                    kbch = 43040;
                    break;
                case C3_4:
                    kbch = 48408;
                    break;
                case C4_5:
                    kbch = 51648;
                    break;
                case C5_6:
                    kbch = 53840;
                    break;
                case C8_9:
                    kbch = 57472;
                    break;
                case C9_10:
                    kbch = 58192;
                    break;                
                default:
                    kbch = 0;
                    break;
            }
        }
        else
        {
            switch (rate)
            {
                case C1_4:
                    kbch = 3072;
                    break;
                case C1_3:
                    kbch = 5232;
                    break;
                case C2_5:
                    kbch = 6312;
                    break;
                case C1_2:
                    kbch = 7032;
                    break;
                case C3_5:
                    kbch = 9552;
                    break;
                case C2_3:
                    kbch = 10632;
                    break;
                case C3_4:
                    kbch = 11712;
                    break;
                case C4_5:
                    kbch = 12432;
                    break;
                case C5_6:
                    kbch = 13152;
                    break;
                case C8_9:
                    kbch = 14232;
                    break;
                case C9_10:
                    fprintf(stderr, "9/10 code rate not supported for DVB-S2 short FECFRAME.\n");
                    exit(1);
                    break;                
                default:
                    kbch = 0;
                    break;
            }
        }
        init_bb_randomiser();
        //set_output_multiple(kbch);
    }

    /*
     * Our virtual destructor.
     */
    bbscrambler_bb_impl::~bbscrambler_bb_impl()
    {
    }

void bbscrambler_bb_impl::init_bb_randomiser(void)
{
    int sr = 0x4A80;
    for (int i = 0; i < FRAME_SIZE_NORMAL; i++)
    {
        int b = ((sr) ^ (sr >> 1)) & 1;
        bb_randomise[i] = b;
        sr >>= 1;
        if(b) sr |= 0x4000;
    }
}

    int
    bbscrambler_bb_impl::work(int noutput_items,
			  const void *input_items,
			  void *output_items)
    {
        const unsigned char *in = (const unsigned char *) input_items;
        unsigned char *out = (unsigned char *) output_items;

        for (int i = 0; i < noutput_items; i += kbch)
        {
            for (int j = 0; j < (int)kbch; ++j)
            {
                out[i + j] = in[i + j] ^ bb_randomise[j];
            }
        }

        // Tell runtime system how many output items we produced.
        return noutput_items;
    }
