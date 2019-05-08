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

#ifdef HAVE_CONFIG_H 
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "demod_sync_cc_impl.h"
#include <gnuradio/expj.h>
#include <math.h>
#include <string.h>
#include <volk/volk.h>

namespace gr {
  namespace fft {

    demod_sync_cc::sptr
    demod_sync_cc::make(int overrate,int fft_size,bool forward,const std::vector<float> &window,bool shift,int nthreads)
    {
      return gnuradio::get_initial_sptr
        (new demod_sync_cc_impl(overrate, fft_size, forward, window, shift, nthreads));
    }

    /*
     * The private constructor
     */
    demod_sync_cc_impl::demod_sync_cc_impl(int overrate,int fft_size,bool forward,const std::vector<float> &window,bool shift,int nthreads)
      : gr::sync_block("demod_sync_cc",
              gr::io_signature::make(1, 1, overrate*fft_size * sizeof(gr_complex)),
              gr::io_signature::make(1, 1, overrate*fft_size * sizeof(gr_complex))),
      d_fft_size(fft_size), d_forward(forward), d_shift(shift),d_overrate(overrate)
    {
	d_fft = new fft_complex(d_fft_size, forward, nthreads);
      if(!set_window(window))
        throw std::runtime_error("demod_cc: window not the same length as fft_size\n");
      if(d_overrate<0)
	throw std::runtime_error("demod_cc: overrate should be greater than 0\n");
	}

    /*
     * Our virtual destructor.
     */
    demod_sync_cc_impl::~demod_sync_cc_impl()
    {
	delete d_fft;
    }
 void
    demod_sync_cc_impl::set_nthreads(int n)
    {
      d_fft->set_nthreads(n);
    }

    int
    demod_sync_cc_impl::nthreads() const
    {
      return d_fft->nthreads();
    }

    bool
    demod_sync_cc_impl::set_window(const std::vector<float> &window)
    {
      if(window.size()==0 || window.size()==d_fft_size) {
    d_window=window;
    return true;
      }
      else
    return false;
    }

    int
    demod_sync_cc_impl::work(int noutput_items,
        gr_vector_const_void_star &input_items,
        gr_vector_void_star &output_items)
    {
      const gr_complex *in = (const gr_complex *) input_items[0];
      gr_complex *out = (gr_complex *) output_items[0];
	
      int noutput_data = noutput_items*d_fft_size;
      gr_complex *in_temp = new gr_complex[noutput_data*d_overrate];
      gr_complex *out_temp = new gr_complex[noutput_data*d_overrate];
      gr_complex *in_temp_ptr = in_temp;
      gr_complex *out_temp_ptr = out_temp;
      if(in_temp!=NULL&&out_temp!=NULL){
	//GR_LOG_WARN(d_logger, "The gr::digital::new rand.");
	//std::cout<<"items:"<<in_temp<<"!";
	//throw std::runtime_error("allocata error");
	}else{
	throw std::runtime_error("allocata error");
	}
      
     for(int i=0;i<d_overrate;i++){
	for(int j=0;j<noutput_data;j++){
	   in_temp[i*noutput_data+j] = in[j*d_overrate+i]; 	
	}	
	}
	/*for(int j=0;j<noutput_data*d_overrate;j++){
	   in_temp[j] = in[j]; 	
	}*/

	unsigned int input_data_size = input_signature()->sizeof_stream_item (0)/d_overrate;
      unsigned int output_data_size = output_signature()->sizeof_stream_item (0);

      int count = 0;
     // gr_complex sum_temp = 0;

      while(count++ < noutput_items*d_overrate) {

      // copy input into optimally aligned buffer
      if(d_window.size()) {
        gr_complex *dst = d_fft->get_inbuf();
        if(!d_forward && d_shift) { 
          unsigned int offset = (!d_forward && d_shift)?(d_fft_size/2):0;
          int fft_m_offset = d_fft_size - offset;
          volk_32fc_32f_multiply_32fc(&dst[fft_m_offset], &in_temp_ptr[0], &d_window[0], offset);
          volk_32fc_32f_multiply_32fc(&dst[0], &in_temp_ptr[offset], &d_window[offset], d_fft_size-offset);
        }
        else {
          volk_32fc_32f_multiply_32fc(&dst[0], in_temp_ptr, &d_window[0], d_fft_size);
        }
      }
      else {
        if(!d_forward && d_shift) {  // apply an ifft shift on the data
          gr_complex *dst = d_fft->get_inbuf();
          unsigned int len = (unsigned int)(floor(d_fft_size/2.0)); // half length of complex array
          memcpy(&dst[0], &in_temp_ptr[len], sizeof(gr_complex)*(d_fft_size - len));
          memcpy(&dst[d_fft_size - len], &in_temp_ptr[0], sizeof(gr_complex)*len);
        }
        else {
          memcpy(d_fft->get_inbuf(), in_temp_ptr, input_data_size);
        }
      }

      // compute the fft
      d_fft->execute();

      // copy result to our output
      if(d_forward && d_shift) {  // apply a fft shift on the data
        unsigned int len = (unsigned int)(ceil(d_fft_size/2.0));
        memcpy(&out_temp_ptr[0], &d_fft->get_outbuf()[len], sizeof(gr_complex)*(d_fft_size - len));
        memcpy(&out_temp_ptr[d_fft_size - len], &d_fft->get_outbuf()[0], sizeof(gr_complex)*len);
      }
      else {
        memcpy (out_temp_ptr, d_fft->get_outbuf (), output_data_size);
      }

      in_temp_ptr  += d_fft_size;
      out_temp_ptr += d_fft_size;
      }
	
      int over_fft_length = d_fft_size*d_overrate;
	  for(int i=0;i<noutput_items;i++){
	  	for(unsigned int j=0;j<d_fft_size;j++){
	  		for(int k=0;k<d_overrate;k++){
	  			out[i*over_fft_length+k*d_fft_size+j] = out_temp[k*noutput_data+i*d_fft_size+j]*gr_expj(-(2*M_PI*(int)(j-d_fft_size/2)*k)/d_fft_size/d_overrate);
	  		}
	  	}
	  }
      // Do <+signal processing+>
      // Tell runtime system how many input items we consumed on
      // each input stream.
	delete []in_temp;
	delete []out_temp;


      // Tell runtime system how many output items we produced.
      return noutput_items;
    }

  } /* namespace fft */
} /* namespace gr */

