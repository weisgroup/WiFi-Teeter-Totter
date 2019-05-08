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


#ifndef INCLUDED_DIGITAL_OFDM_FRAME_OS_ACQUISITION_H
#define INCLUDED_DIGITAL_OFDM_FRAME_OS_ACQUISITION_H

#include <gnuradio/digital/api.h>
#include <gnuradio/sync_block.h>
#include <vector>

namespace gr {
  namespace digital {

    /*!
     * \brief <+description of block+>
     * \ingroup digital
     *
     */
    class DIGITAL_API ofdm_frame_os_acquisition : virtual public gr::sync_block
    {
     public:
      typedef boost::shared_ptr<ofdm_frame_os_acquisition> sptr;

      /*! 
       * Make an OFDM correlator and equalizer.
       *
       * \param occupied_carriers   The number of subcarriers with data in the received symbol
       * \param fft_length          The size of the FFT vector (occupied_carriers + unused carriers)
       * \param cplen		    The length of the cycle prefix
       * \param known_symbol        A vector of complex numbers representing a known symbol at the
       *                            start of a frame (usually a BPSK PN sequence)
       * \param max_fft_shift_len   Set's the maximum distance you can look between bins for correlation
       */
      static sptr make(unsigned int overrate,unsigned int occupied_carriers, unsigned int fft_length,unsigned int cplen,const std::vector<gr_complex> &known_symbol, unsigned int max_fft_shift_len=4);
      /*!
       * \brief Return an estimate of the SNR of the channel
       */
      virtual float snr() = 0;
    };

  } // namespace digital
} // namespace gr

#endif /* INCLUDED_DIGITAL_OFDM_FRAME_OS_ACQUISITION_H */

