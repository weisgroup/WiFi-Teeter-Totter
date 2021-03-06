#!/usr/bin/env python
#
# Copyright 2006,2007,2011,2013 Free Software Foundation, Inc.
# 
# This file is part of GNU Radio
# 
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# 
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
# 

from gnuradio import gr
from gnuradio import eng_notation
from gnuradio.eng_option import eng_option
from optparse import OptionParser

from gnuradio import blocks
from gnuradio import digital
from gnuradio import filter

# from current dir
from receive_path import receive_path
from uhd_interface import uhd_receiver

import struct, sys

class my_top_block(gr.top_block):
    def __init__(self, callback, options):
        gr.top_block.__init__(self)

        if(options.rx_freq is not None):
            self.source = uhd_receiver(options.args,
                                       options.bandwidth, options.rx_freq, 
                                       options.lo_offset, options.rx_gain,
                                       options.spec, options.antenna,
                                       options.clock_source, options.verbose)
        elif(options.from_file is not None):
            self.source = blocks.file_source(gr.sizeof_gr_complex, options.from_file)
        else:
            self.source = blocks.null_source(gr.sizeof_gr_complex)

        # Set up receive path
        # do this after for any adjustments to the options that may
        # occur in the sinks (specifically the UHD sink)

        self.rxpath = receive_path(callback, options)

        self.connect(self.source, self.rxpath)

        # self.resampler = filter.fractional_resampler_cc(0, options.overrate)
        # options.overrate = 1
        # self.rxpath = receive_path(callback, options)

        # self.connect(self.source, self.resampler, self.rxpath)

# /////////////////////////////////////////////////////////////////////////////
#                                   main
# /////////////////////////////////////////////////////////////////////////////

def main():

    global n_rcvd, n_right,n_cbits,n_tbits, n_bytes, n_right_bytes
        
    n_rcvd = 0
    n_right = 0
    n_cbits = 0
    n_tbits = 0
    n_bytes = 0
    n_right_bytes = 0
    if 1:
        file_handler = open('rx_payload.txt', 'w')
        file_handler.write('')
        file_handler.close()

    def rx_callback(ok, payload,cbits,tbits):
        global n_rcvd, n_right,n_cbits,n_tbits, n_bytes, n_right_bytes

        n_rcvd += 1
        n_tbits+=tbits
        n_cbits+=cbits
        ber = float(100.0*float(n_tbits-n_cbits)/float(n_tbits))
        ber_p = float(100.0*float(tbits-cbits)/float(tbits))
    
        (pktno,) = struct.unpack('!H', payload[0:2])

        if ok:
            n_right += 1
        print "ok: %r \t pktno: %d \t n_rcvd: %d \t n_right: %d\t ber:%.4f%%\t "% (ok, pktno, n_rcvd, n_right,ber)
        # print "ok: %r \t pktno: %d \t n_rcvd: %d \t n_right: %d" % (ok, pktno, n_rcvd, n_right)
        if 1:
            printlst = list()
            byte_payload_num = 0
            byte_right_num = 0
            printlst.append(hex(pktno))
            for x in payload[2:]:
                t = hex(ord(x)).replace('0x', '')
                
                #if t == hex(pktno & 0xff).replace('0x', ''):
                if t == '0':
                    byte_right_num += 1
                byte_payload_num += 1
                if(len(t) == 1):
                    t = '0' + t
                printlst.append(t)
            printable = ''.join(printlst)

            n_bytes += byte_payload_num 
            n_right_bytes += byte_right_num

            # write the payload to a file
            if 0:
                file_handler = open('rx_payload.txt', 'a')
                file_handler.write(printable +'\n')
                file_handler.close()

            #print "pktno: %d \t right byte: %d \t byte accep rate: %.4f%%" % ( pktno, byte_right_num, float(n_right_bytes)/n_bytes * 100)

    parser = OptionParser(option_class=eng_option, conflict_handler="resolve")
    expert_grp = parser.add_option_group("Expert")
    parser.add_option("","--discontinuous", action="store_true", default=False,
                      help="enable discontinuous")
    parser.add_option("","--from-file", default=None,
                      help="input file of samples to demod")

    receive_path.add_options(parser, expert_grp)
    uhd_receiver.add_options(parser)
    digital.ofdm_demod.add_options(parser, expert_grp)

    (options, args) = parser.parse_args ()

    if options.from_file is None:
        if options.rx_freq is None:
            sys.stderr.write("You must specify -f FREQ or --freq FREQ\n")
            parser.print_help(sys.stderr)
            sys.exit(1)

    # build the graph
    tb = my_top_block(rx_callback, options)

    r = gr.enable_realtime_scheduling()
    if r != gr.RT_OK:
        print "Warning: failed to enable realtime scheduling"

    tb.start()                      # start flow graph
    tb.wait()                       # wait for it to finish

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass
