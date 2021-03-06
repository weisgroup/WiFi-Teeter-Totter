/* -*- c++ -*- */
/* 
 * Copyright 2017 Free Software Foundation, Inc.
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "ofdm_oversampler_impl.h"
#include <gnuradio/expj.h>
#include <cstdio>


namespace gr {
  namespace digital {

    ofdm_oversampler::sptr
    ofdm_oversampler::make(int overrate,unsigned int fft_length,unsigned int symbol_length,unsigned int timeout)
    {
      return gnuradio::get_initial_sptr
        (new ofdm_oversampler_impl(overrate, fft_length, symbol_length, timeout));
    }

    /*
     * The private constructor
     */
    ofdm_oversampler_impl::ofdm_oversampler_impl(int overrate,unsigned int fft_length,unsigned int symbol_length,unsigned int timeout)
      : gr::block("ofdm_oversampler",
              gr::io_signature::make2(2, 2, sizeof(gr_complex), sizeof(char)),
              gr::io_signature::make2(2, 2, sizeof(gr_complex)*fft_length*overrate, sizeof(char)*fft_length)),
	d_state(STATE_NO_SIG), d_timeout_max(timeout),
	d_fft_length(fft_length), d_symbol_length(symbol_length),d_overrate(overrate)
    {
		GR_LOG_WARN(d_logger, "The gr::digital::ofdm_oversampler block has been deprecated.");

      set_relative_rate(1.0/(double) fft_length);   // buffer allocator hint
	}

    /*
     * Our virtual destructor.
     */
    ofdm_oversampler_impl::~ofdm_oversampler_impl()
    {
    }

    void
    ofdm_oversampler_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      // FIXME do we need more
      //int nreqd  = (noutput_items-1) * d_symbol_length + d_fft_length;
      int nreqd  = (d_symbol_length + d_fft_length)*d_overrate;
      unsigned ninputs = ninput_items_required.size();
      for(unsigned i = 0; i < ninputs; i++)
	ninput_items_required[i] = nreqd;
    }

    int
    ofdm_oversampler_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
	//std::cout<<"im fine";
      const gr_complex *iptr = (const gr_complex*)input_items[0];
      const char *trigger = (const char*)input_items[1];

      gr_complex *optr = (gr_complex*)output_items[0];
      char *outsig = (char*)output_items[1];

      //FIXME: we only process a single OFDM symbol at a time; after the preamble, we can
      // process a few at a time as long as we always look out for the next preamble.

      unsigned int index = d_fft_length*d_overrate+4;  // start one fft length into the input so we can always look back this far

      outsig[0] = 0; // set output to no signal by default

      // Search for a preamble trigger signal during the next symbol length
      while((d_state != STATE_PREAMBLE) && (index <= (d_symbol_length+d_fft_length)*d_overrate)) {
	if(trigger[index]) {
	  //index = index-1;
	  outsig[0] = 1; // tell the next block there is a preamble coming
	  d_state = STATE_PREAMBLE;
	}
	else
	  index++;
      }

      unsigned int i, pos, ret;
      switch(d_state) {
      case(STATE_PREAMBLE):
	// When we found a preamble trigger, get it and set the symbol boundary here
	for(i = (index - d_fft_length*d_overrate + 1); i <= index; i++) {
	  *optr++ = iptr[i-3];
	}

	d_timeout = d_timeout_max; // tell the system to expect at least this many symbols for a frame
	d_state = STATE_FRAME;
	consume_each(index - d_fft_length*d_overrate + 1); // consume up to one fft_length away to keep the history
	ret = 1;
	break;

      case(STATE_FRAME):
	// use this state when we have processed a preamble and are getting the rest of the frames
	//FIXME: we could also have a power squelch system here to enter STATE_NO_SIG if no power is received

	// skip over fft length history and cyclic prefix
	pos = d_symbol_length*d_overrate;         // keeps track of where we are in the input buffer
	while(pos < (d_symbol_length + d_fft_length)*d_overrate) {
	  *optr++ = iptr[pos++ -3];
	}

	if(d_timeout-- == 0) {
	  printf("TIMEOUT\n");
	  d_state = STATE_NO_SIG;
	}

	consume_each(d_symbol_length*d_overrate); // jump up by 1 fft length and the cyclic prefix length
	ret = 1;
	break;

      case(STATE_NO_SIG):
      default:
	consume_each(index-d_fft_length*d_overrate); // consume everything we've gone through so far leaving the fft length history
      ret = 0;
      break;
      }

      return ret;
    }

  } /* namespace digital */
} /* namespace gr */

