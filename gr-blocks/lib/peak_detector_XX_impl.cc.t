/* -*- c++ -*- */
/*
 * Copyright 2007,2010,2013 Free Software Foundation, Inc.
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

// @WARNING@

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "@NAME_IMPL@.h"
#include <gnuradio/io_signature.h>
#include <string.h>
#include <limits>

namespace gr {
  namespace blocks {

    @NAME@::sptr
    @NAME@::make(float threshold_factor_rise,
                 float threshold_factor_fall,
                 int look_ahead, float alpha)
    {
      return gnuradio::get_initial_sptr
        (new @NAME_IMPL@(threshold_factor_rise,
                         threshold_factor_fall,
                         look_ahead, alpha));
    }

    @NAME_IMPL@::@NAME_IMPL@(float threshold_factor_rise,
                             float threshold_factor_fall,
                             int look_ahead, float alpha)
    : sync_block("@BASE_NAME@",
                    io_signature::make(1, 1, sizeof(@I_TYPE@)),
                    io_signature::make(1, 1, sizeof(char))),
      d_fft_length(look_ahead), d_threshold(alpha),d_peak_val(0),d_not_found(0),d_times_ref(50)
    {
    }

    @NAME_IMPL@::~@NAME_IMPL@()
    {
    }

    int
    @NAME_IMPL@::work(int noutput_items,
                      gr_vector_const_void_star &input_items,
                      gr_vector_void_star &output_items)
    {

      @I_TYPE@ *iptr = (@I_TYPE@*)input_items[0];
      char *optr = (char *)output_items[0];
      
      memset(optr, 0, noutput_items*sizeof(char));
      int i = 0;
      float total_val = 0;
      float ave_val = 0;
      float peak_val=0;
      int peak_index=0;
      float times_ref=0;
	//std::cout<<"items:"<<noutput_items<<"!";

      for(i=0;i<noutput_items;i++){
      	total_val = total_val + iptr[i];
	if(peak_val<iptr[i]){
      		peak_val=iptr[i];
      	}
      }
      	
      ave_val = total_val/noutput_items;
      times_ref =peak_val/ave_val;
      if(times_ref>d_times_ref&&peak_val>(0.5-d_not_found*0.35)*d_peak_val&&(ave_val>0.1*d_threshold)){
	//d_times_ref=0.7*times_ref;
	//if(d_times_ref/times_ref>0.7||d_times_ref/times_ref<0.2){
	//d_times_ref=d_times_ref+(0.7*times_ref-d_times_ref)*0.1;}
	d_not_found=0;
      d_peak_val = peak_val; 
      i = 0;

      while(i<noutput_items){
      	if((iptr[i]/ave_val)>(0.5*times_ref)&&(iptr[i]>d_threshold)){
      		peak_val = iptr[i];
      		peak_index = i;
      		if(i+(d_fft_length/2)+1<noutput_items){
      			for(int j=0;j<=(d_fft_length/2);j++){
      				if(iptr[i+j]>peak_val){
      					peak_val = iptr[i+j];
      					peak_index = i+j;
      				}
      			}
		optr[peak_index-7]=1;
      		i=i+(d_fft_length/2);
      		}else{
		for(int j=0;j<=(noutput_items-i-2);j++){
      				if(iptr[i+j]>peak_val){
      					peak_val = iptr[i+j];
      					peak_index = i+j;
      				}
      			}
		optr[peak_index-7]=1;
      		i=noutput_items;
		}
      	}else{
		i=i+1;
	}
      }}else{
	d_not_found+=noutput_items/(d_fft_length*19*5/4);
	}
      return noutput_items;
    }

  } /* namespace blocks */
} /* namespace gr */
