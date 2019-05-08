### Wi-Fi Teeter-Totter

Wi-Fi Teeter-Totter PHY layer is built atop open source GNURadio <http://gnuradio.org/>

### Introduction

This is the source code of Overclocking [1]. The source code is maintained locally by our research group and not published until the related paper has been accepted recently.

The main idea of overclocking is that the OFDM receiver takes full power of its high clock rate to boost the decoding ability so that decode the signal in low SNR.

Overclocking is implemented in PHY layer, which based on the standard of IEEE 802.11n. It is realized by the revising of OFDM example “benchmark_rx” in GNURadio/USRP B210. 

GNU Radio is an open source development toolkit that provides the digital signal processing (DSP) blocks to implement software defined radios. To ensure the runtime performance, the DSP blocks are implemented using C++. Through exploiting Python to connect different DSP blocks, we can fast emulate and prototype the low-level functionality of wireless protocol families (e.g. Wi-Fi, RFID, Bluetooth and Cellular, etc.).

### Environment and Installation

* Our project is build on the Ubuntu 16.04 with GNURadio 3.7.11, and you should have installed the dependent package and USRP hardware driver beforehand.

* Installed GNUradio 3.7.11 in your Linux PC 

  Related guide: <http://wikignuradio.org/index.php/InstallingGR>

* Compile the code we uploaded using the commands below:

```shell
$ cd ../gnuradio-3.7.11/build (go to the build file of gnuradio)
$ cmake ../
$make
$sudo make install
```

After finishing these steps,you can test it in your USRP platform.

### Usage

Transmitter: 

```shell
./benchmark_tx.py -f 930e6 -fft-length 64 -occupied-tones 48 –cp-length 16 –log –v –m bpsk –bandwidth 500000
```

Receiver: 

```shell
 ./benchmark_rx.py -f 930e6 -fft-length 64 -occupied-tones 48 –cp-length 16 –log –v –m bpsk -A TX/RX –bandwidth 1000000 -overrate 2
```

The commands above serves transmitter and receiver respectively, and overrate represent the overclocking rate in the communication system. 

### Reference

[1]*”Wi-Fi Teeter-Totter: Overclocking OFDM for Internet of Things”,IEEE INFOCOM(2018)*