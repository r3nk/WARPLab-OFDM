WARPLab-OFDM
============

This package contains an [OFDM][] framework that facilitates physical transmissions by means of [WARP][], the Wireless Open-Access Research Platform, and [WARPLab][], a framework that allows to access WARP hardware from the [MATLAB][] workspace. 
The framework may serve as a basis for prototyping [PHY][] mechanisms based on OFDM. 
Besides physical transmissions with WARP, the framework also allows for a simulated [AWGN][] channel by MATLAB and may easily be extended by custom wireless channel models and interfaces of other [SDR][] platforms. 



Usage
------------

OFDM transmissions are initiated by the MATLAB script `ofdm_run`. 
It allows to set up all PHY parameters and switches for visualization, such as plots and scatter diagrams. 
Furthermore, it allows for execution in a loop and to save recorded transmissions to files. 

The `ofdm_run` script calls the `ofdm_func()` function, which implements the actual OFDM functionality. 
The latter extracts the passed parameters, generates the OFDM signal, passes it to the transceiver for transmission, demodulates the received signals, and calculates error rates. 
The generated OFDM signal consists of a sequence of `Npil` pilot symbols, followed by a sequence of `Ndat` data symbols. 
Pilot symbols are employed to estimate the channel coefficients as well as the corresponding SNRs. 
Furthermore, pilot symbols are used to compensate the Carrier Frequency Offset (CFO) in practical setups, such as WARPLab, prior to channel estimation. 
The channel used to transmit the  generated OFDM signal can be configured by the parameter `params.transceiver` in `ofdm_run`, which also allows to set up transmissions by means of WARPLab. 

The functions `warplab_2x2_multi_transceive()` and `warplab_4x4_multi_transceive()` transmit a block of samples via WARPLab from a first WARP node to a second WARP node, where the received signals are captured. 
These functions additionally perform preamble synchronization in order to return the vectors of received samples in sync with the transmitted samples. 
The vectors of received samples thereby exactly match the transmitted vectors in length, which facilitates to abstract from the physical channel by the application. 
Note that these functions are not tied to OFDM but may be used to transceive any kinds of signals between WARPs. 



License
------------

Please see [LICENSE.txt](https://github.com/r3nk/WARPLab-OFDM/blob/master/LICENSE.txt) for detailed licensing information. 
Note that any use of the framework which results in an academic publication or other publication which includes a bibliography must include a citation to the author's publication *A Rapid Prototyping Framework for Practical OFDMA Systems using Software Defined Radios*. 
A BibTeX cite key is available at the [SEEMOO][SEEMOO-OFDM] page. 



Compatibility
------------

The framework was developed and tested with the following configuration:

*   WARP FPGA v2.2
*   WARPLab v5.02
*   MATLAB R2007b



About 
------------

The OFDM framework was devised at the [Secure Mobile Networking Lab][SEEMOO] at [Technische Universität Darmstadt][TUD] in spring 2012 as part of the author's Master thesis with the title *Dynamic Subchannel Allocation in OFDMA-Based Wireless Mesh Networks*. 
The goal of the thesis was to quantify multi-user channel diversity gains by assigning OFDMA subchannels dynamically between communication links of several concurrent transmitters and receivers as a function of measured channel conditions in a decentralized topology.  
The OFDM framework provided in this package is a stripped down variant of the more complex and more specific OFDMA framework devised in the thesis. 
It is meant to serve as a generic foundation for new innovative designs based on OFDM. 



References
------------

This OFDM framework at GitHub:  
<https://github.com/r3nk/WARPLab-OFDM>

Secure Mobile Networking Lab - SEEMOO:  
<http://www.seemoo.tu-darmstadt.de>

Rice University WARP Project:  
<http://warp.rice.edu>


[OFDM]: http://en.wikipedia.org/wiki/OFDM "Orthogonal Frequency-Division Multiplexing"
[WARP]: http://warp.rice.edu "Rice University WARP Project"
[WARPLab]: http://warp.rice.edu/trac/wiki/WARPLab "WARPLab Framework"
[MATLAB]: http://www.mathworks.com/products/matlab "MathWorks MATLAB"
[PHY]: http://en.wikipedia.org/wiki/PHY "Physical Layer"
[AWGN]: http://en.wikipedia.org/wiki/AWGN "Additive White Gaussian Noise"
[SDR]: http://en.wikipedia.org/wiki/Software-defined_radio "Software-Defined Radio"
[SEEMOO-OFDM]: http://www.seemoo.tu-darmstadt.de/ofdm "OFDM framework at SEEMOO"
[SEEMOO]: http://www.seemoo.tu-darmstadt.de "Secure Mobile Networking Lab"
[TUD]: http://www.tu-darmstadt.de "Technische Universität Darmstadt"
