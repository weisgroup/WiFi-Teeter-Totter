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


#ifndef INCLUDED_FFT_DEMOD_SYNC_CC_H
#define INCLUDED_FFT_DEMOD_SYNC_CC_H

#include <gnuradio/fft/api.h>
#include <gnuradio/sync_block.h>

namespace gr {
  namespace fft {

    /*!
     * \brief <+description of block+>
     * \ingroup fft
     *
     */
    class FFT_API demod_sync_cc : virtual public gr::sync_block
    {
     public:
      typedef boost::shared_ptr<demod_sync_cc> sptr;

      virtual void set_nthreads(int n) = 0;

      virtual int nthreads() const = 0;

      virtual bool set_window(const std::vector<float> &window) = 0;
      static sptr make(int overrate,int fft_size,bool forward,const std::vector<float> &window,bool shift=false,int nthreads=1);
    };

  } // namespace fft
} // namespace gr

#endif /* INCLUDED_FFT_DEMOD_SYNC_CC_H */

