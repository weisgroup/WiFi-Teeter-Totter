/* -*- c++ -*- */
/* 
 * Copyright 2018 Free Software Foundation, Inc.
 * 
 * This file is part of GNU Radio
 * 
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_FFT_DEMOD_SYNC_CC_IMPL_H
#define INCLUDED_FFT_DEMOD_SYNC_CC_IMPL_H

#include <gnuradio/fft/demod_sync_cc.h>
#include <gnuradio/fft/fft.h>

namespace gr {
  namespace fft {

    class demod_sync_cc_impl : public demod_sync_cc
    {
     private:
	fft_complex          *d_fft;
	//gr_complex		*in_temp;
	//gr_complex		*out_temp;
      	unsigned int          d_fft_size;
      	std::vector<float>    d_window;
      	bool                  d_forward;
      	bool                  d_shift;
	int		      d_overrate;

     public:
      demod_sync_cc_impl(int overrate,int fft_size,bool forward,const std::vector<float> &window,bool shift,int nthreads);
      ~demod_sync_cc_impl();

      void set_nthreads(int n);
      int nthreads() const;
      bool set_window(const std::vector<float> &window);

      int work(int noutput_items,
         gr_vector_const_void_star &input_items,
         gr_vector_void_star &output_items);
    };

  } // namespace fft
} // namespace gr

#endif /* INCLUDED_FFT_DEMOD_SYNC_CC_IMPL_H */

