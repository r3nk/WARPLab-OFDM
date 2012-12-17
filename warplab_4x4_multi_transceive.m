%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OFDM Framework
%
% Copyright (c) 2012 Robin Klose 
%
% Parts of the concepts regarding preamble synchronization are taken from
% the WARPLab workshop exercises available at: http://warp.rice.edu
%
% WARP and WARPLab are licensed under the 
% RICE UNIVERSITY SOFTWARE DISTRIBUTION LICENSE
% 
% [Wireless Open-Access Research Platform Designs & Implementations] 
% (the "Software")
% 
% Copyright (c) 2006-2012, Rice University.  All rights reserved.
% 
% The Software is available for download and use subject to
% the terms and conditions of this License.  Access or use of
% the Software constitutes acceptance and agreement to the
% terms and conditions of this License.
% 
% Redistribution and use of the Software in source and binary forms,
% with or without modification, are permitted provided that the
% following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright
% notice, this list of conditions and the capitalized paragraph below.
% 
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the capitalized paragraph below
% in the documentation and/or other materials provided with the
% distribution. 
% 
% 3. Any use of The Software or the underlying hardware platform which
% results in an academic publication or other publication which includes a
% bibliography must include a citation to the WARP project:
% "Rice University WARP Project, http://warp.rice.edu"
% 
% 4. Except as required to comply with condition #3 above, the names of
% Rice University or its faculty, staff or students may not be used to
% endorse or promote products derived from the Software without specific
% prior written permission.
% 
% THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS, IMPLIED OR
% STATUTORY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, WARRANTIES OF
% ACCURACY, COMPLETENESS, NONINFRINGEMENT, MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED.  ACCESS OR USE OF THE SOFTWARE IS
% ENTIRELY AT THE USER'S RISK.  IN NO EVENT SHALL RICE UNIVERSITY OR ITS
% FACULTY, STAFF OR STUDENTS BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
% NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
% THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
% THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  THE
% SOFTWARE USER SHALL INDEMNIFY, DEFEND AND HOLD HARMLESS RICE UNIVERSITY
% AND ITS FACULTY, STAFF AND STUDENTS FROM ANY AND ALL CLAIMS, ACTIONS,
% DAMAGES, LOSSES, LIABILITIES, COSTS AND EXPENSES, INCLUDING ATTORNEYS'
% FEES AND COURT COSTS, DIRECTLY OR INDIRECTLY ARISING OUR OF OR IN
% CONNECTION WITH ACCESS OR USE OF THE SOFTWARE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit a block of samples via WARPLab 4x4 and receive it at a second
% node. This code allows for simultaneous transmissions with up to four
% radios on WARP board 1 and to receive the transmitted data with up to
% four radios on WARP board 2. The WARP boards should be running the
% WARPLab 4x4 reference design. 
% TX1 is node 1 radio 1. TX2 is node 1 radio 2. 
% TX3 is node 1 radio 3. TX4 is node 1 radio 4. 
% RX1 is node 2 radio 1. RX2 is node 2 radio 2. 
% RX3 is node 2 radio 3. RX4 is node 2 radio 4. 
%
% Provide a structure with the following parameters: 
%
% n_tx:        Number of transmit antennas. Should be in range [1,4]. 
% n_rx:        Number of receive antennas.  Should be in range [1,4]. 
% tx_block:    Data to transmit. Should consist of n_tx column vectors.
%              If it consists of more column vectors, these additional
%              vectors will be transmitted as is without oversampling and
%              without being shifted to an intermediate frequency. 
%              If oversampling is applied to the transmission data, these
%              additional vectors will be repeated 'oversamp' times
%              instead. This feature is designed to be used for artificial
%              noise interferers. 
% oversamp:    An oversampling pulse shaping filter is applied on the data
%              for transmission if the oversamp parameter is greater 1. 
%              If this parameter is 0 or 1, the data for transmission is
%              not processed by the pulse shaping filter and transmitted as
%              is.
% int_f:       Intermediate frequency in the baseband. May be used to avoid
%              DC attenuation. Keep in mind that the baseband signal's
%              bandwidth should not exceed 9.5 MHz in nominal mode. 
%              For and intermediate frequency of 5 MHz, for instance, pass:
%              int_f = 5e6; 
% channel:     [1,14] The carrier channel to be used. 
%              Pass 0 to not change a previously setup channel.
% tx1_gain_bb: [0, 3] WARPLab TX1 baseband gain.
% tx2_gain_bb: [0, 3] WARPLab TX2 baseband gain.
% tx3_gain_bb: [0, 3] WARPLab TX3 baseband gain.
% tx4_gain_bb: [0, 3] WARPLab TX4 baseband gain.
% tx1_gain_rf: [0,63] WARPLab TX1 radio frequency gain.
% tx2_gain_rf: [0,63] WARPLab TX2 radio frequency gain.
% tx3_gain_rf: [0,63] WARPLab TX3 radio frequency gain.
% tx4_gain_rf: [0,63] WARPLab TX4 radio frequency gain.
% rx1_gain_bb: [0,31] WARPLab RX1 baseband gain.
% rx2_gain_bb: [0,31] WARPLab RX2 baseband gain.
% rx3_gain_bb: [0,31] WARPLab RX3 baseband gain.
% rx4_gain_bb: [0,31] WARPLab RX4 baseband gain.
% rx1_gain_rf: [1, 3] WARPLab RX1 radio frequency gain.
% rx2_gain_rf: [1, 3] WARPLab RX2 radio frequency gain.
% rx3_gain_rf: [1, 3] WARPLab RX3 radio frequency gain.
% rx4_gain_rf: [1, 3] WARPLab RX4 radio frequency gain.
% samp_offset: Manually fine tune the alignment of the block of received
%              samples. 
% print_debug: Print debug messages. 
% plot_filter: Switch to plot impulse response of SRRC filter. 
% plot_spectrum: Switch to plot the spectrum(s) of the received signal(s). 
% plot_waveform: Switch to plot the waveform(s) of the received signal(s). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rx_block = warplab_4x4_multi_transceive(params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_tx        = params.n_tx;
n_rx        = params.n_rx;
tx_block    = params.tx_block;
oversamp    = params.oversamp;
int_f       = params.int_f;
channel     = params.channel;
tx1_gain_bb = params.tx1_gain_bb;
tx2_gain_bb = params.tx2_gain_bb;
tx3_gain_bb = params.tx3_gain_bb;
tx4_gain_bb = params.tx4_gain_bb;
tx1_gain_rf = params.tx1_gain_rf;
tx2_gain_rf = params.tx2_gain_rf;
tx3_gain_rf = params.tx3_gain_rf;
tx4_gain_rf = params.tx4_gain_rf;
rx1_gain_bb = params.rx1_gain_bb;
rx2_gain_bb = params.rx2_gain_bb;
rx3_gain_bb = params.rx3_gain_bb;
rx4_gain_bb = params.rx4_gain_bb;
rx1_gain_rf = params.rx1_gain_rf;
rx2_gain_rf = params.rx2_gain_rf;
rx3_gain_rf = params.rx3_gain_rf;
rx4_gain_rf = params.rx4_gain_rf;
samp_offset = params.samp_offset;
print_debug = params.print_debug;
plot_filter = params.plot_filter;
plot_spectrum = params.plot_spectrum;
plot_waveform = params.plot_waveform;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assumption of maximum jitter in the sync trigger as a number of samples:
max_sync_jitter = 100;

% Sampling frequency and time: 
f_sam = 40e6;
t_sam = 1/f_sam;

% Extended number of transmissions (for additional noise tx vectors)
n_tx_ext = size(tx_block, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Debugging
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Print arguments:
if print_debug
    fprintf('n_tx = %d\n', n_tx);
    fprintf('n_rx = %d\n', n_rx);
    fprintf('size(tx_block) = %d, %d\n', ...
            size(tx_block, 1), size(tx_block, 2));
    fprintf('oversamp = %d\n', oversamp);
    fprintf('channel = %d\n', channel);
    fprintf('tx1_gain_bb = %d\n', tx1_gain_bb);
    fprintf('tx2_gain_bb = %d\n', tx2_gain_bb);
    fprintf('tx3_gain_bb = %d\n', tx3_gain_bb);
    fprintf('tx4_gain_bb = %d\n', tx4_gain_bb);
    fprintf('tx1_gain_rf = %d\n', tx1_gain_rf);
    fprintf('tx2_gain_rf = %d\n', tx2_gain_rf);
    fprintf('tx3_gain_rf = %d\n', tx3_gain_rf);
    fprintf('tx4_gain_rf = %d\n', tx4_gain_rf);
    fprintf('rx1_gain_bb = %d\n', rx1_gain_bb);
    fprintf('rx2_gain_bb = %d\n', rx2_gain_bb);
    fprintf('rx3_gain_bb = %d\n', rx3_gain_bb);
    fprintf('rx4_gain_bb = %d\n', rx4_gain_bb);
    fprintf('rx1_gain_rf = %d\n', rx1_gain_rf);
    fprintf('rx2_gain_rf = %d\n', rx2_gain_rf);
    fprintf('rx3_gain_rf = %d\n', rx3_gain_rf);
    fprintf('rx4_gain_rf = %d\n', rx4_gain_rf);
    fprintf('plot_filter = %d\n', plot_filter);
    fprintf('plot_spectrum = %d\n', plot_spectrum);
    fprintf('plot_waveform = %d\n', plot_waveform);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of pulse shaping filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The pulse shaping filter is used to upsample the preamble of
% transmission and thereby limits the bandwidth of the transmitted signal. 
% 
% The preamble will be detected a the receiver to identify the start
% position of the transmission sequence in the received signal. 
% Furthermore, the pulse shaping filter is applied to the data for
% tranmission if the oversamp paramter is greater than 1. 
% The pulse shaping filter is a Square-Root Raised Cosine (SRRC) filter.


% FIR filter order: 
filter_order = 64; 

if oversamp < 1
    oversamp = 1;
end

% Oversampling by filter:
if oversamp > 1
    filter_samp = oversamp; 
else
    filter_samp = 8; 
end

% Group delay: Time between the input to the filter and the filter's peak
% response in terms of samples.
filter_input_delay = filter_order / (filter_samp * 2);
filter_output_delay = filter_order / 2;

% Filter rolloff factor: 
filter_rolloff = 0.3;

% Create the filter: 
rrcfilter = rcosine(1, filter_samp, 'fir/sqrt', ... 
                    filter_rolloff, filter_input_delay); 

% Plot the filter's impulse response in a stem plot: 
if plot_filter
    figure;
    stem(rrcfilter);
    title('Raised Cosine Impulse Response');
    xlabel('Sample'); ylabel('Amplitude');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of preamble and correlation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the preamble. The preamble is a BPSK-modulated Barker sequence: 
preamble = [-1;-1;-1;1;-1;0;0;0;0;0;0;0;0]; 

% The WARPLab transmit buffers can store up to 2^14 samples: 
if (2 * filter_output_delay + ...
    filter_samp * length(preamble) + ...
    oversamp * size(tx_block, 1) + ...
    2 * max_sync_jitter > 2^14)
    
    error('Too large data block!');
end

% Create a reference matrix for detection of the preamble in the received
% signal. The 'corr_window' gives the size of the correlation window as a
% number of samples. The preamble is expected to be found within the
% correlation window. 
% The reference matrix contains in each column a circularly shifted copy 
% of the upsampled preamble. The first column of the reference matrix
% contains the zero-padded upsampled preamble. The i-th column of the
% reference matrix contains the circular shift of the reference samples
% vector. 
% 
preamble_us = upsample(preamble, filter_samp);
preamble_us_len = length(preamble_us);
corr_window = 2 * max_sync_jitter + preamble_us_len;
reference_samples = zeros(corr_window,1); 
reference_samples(1:preamble_us_len) = preamble_us; 
reference_matrix = toeplitz(reference_samples,...
                    circshift(reference_samples(corr_window:-1:1),1)); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concatenate preamble and transmission symbols and upsample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

preamble_filt_len = 2*filter_output_delay + filter_samp*length(preamble);
                % = length(preamble_filt);

% The preamble is inserted in the first transmission vector. If there are
% parallel transmissions, the other vectors will be padded with zeros. 
if oversamp > 1
    % Oversample both the preamble and data: 
    tx_block_pa = [zeros(length(preamble), size(tx_block, 2)); 
                   tx_block];
    tx_block_pa(1:length(preamble), 1) = preamble;
    tx_block_pa_filt = rcosflt(tx_block_pa, 1, oversamp, ... 
                               'filter', rrcfilter); 
    % If more than n_tx vectors are given, treat these additional vectors
    % differently: Instead of oversampling, repeat them 'oversamp' times. 
    % This feature is designed for artificial noise interferers. 
    if (n_tx_ext > n_tx)
        tx_block_pa_filt(preamble_filt_len + 1 : end, n_tx + 1 : end) = ...
            repmat(tx_block(:, n_tx + 1 : end), oversamp, 1);
    end
else
    % Oversample only the preamble: 
    preamble_filt = rcosflt(preamble, 1, filter_samp, 'filter', rrcfilter);
    tx_block_pa_filt = [zeros(preamble_filt_len, size(tx_block, 2));
                        tx_block];
    tx_block_pa_filt(1:preamble_filt_len, 1) = preamble_filt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Upconvert to an intermediate frequency if desired. 
% This may be useful to avoid DC attenuation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if int_f
    time = t_sam * (0 : 1 : size(tx_block_pa_filt, 1) - 1).'; 
    tx_block_pa_filt_up = bsxfun(@times, tx_block_pa_filt, ... 
                                 exp(j*2*pi*int_f*time));
    if (n_tx_ext > n_tx)
        tx_block_pa_filt_up(:, n_tx + 1 : end) = ...
                        tx_block_pa_filt(:, n_tx + 1 : end);
    end
else
    tx_block_pa_filt_up = tx_block_pa_filt;
end

% Normalize the vectors to transmit to [-1, 1] in order to meet the DAC
% requirements of the WARP transmitter. 
tx_block_pa_filt_up_max_i = max(abs(real(tx_block_pa_filt_up(:))));
tx_block_pa_filt_up_max_q = max(abs(imag(tx_block_pa_filt_up(:))));
tx_block_pa_filt_up_max_iq = max(tx_block_pa_filt_up_max_i, ... 
                                 tx_block_pa_filt_up_max_q);

tx_block_pa_filt_up = tx_block_pa_filt_up / tx_block_pa_filt_up_max_iq;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmission via WARPLab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load some global definitions:
warplab_defines

% Create socket handles and initialize nodes:
[socketHandles, packetNum] = warplab_initialize;

% Socket handle for the magic SYNC: 
udp_sync = socketHandles(1);

% Open two handles: Node 1 is used for transmission, node 2 for reception:
udp_node1 = socketHandles(2);
udp_node2 = socketHandles(3);

% Delay the transmission by the maximum assumed jitter in the sync trigger.
% This enables to detect the preamble if the jitter in the sync trigger
% does not exceed 'max_sync_jitter' samples. 
tx_delay = max_sync_jitter;

% Number of samples to transmit: 
tx_length = size(tx_block_pa_filt_up, 1);

% TX mode. 
% 0: singe transmission 
% 1: continuous transmission
tx_mode = 0;

% Download the WARPLab parameters to the WARP nodes: 
warplab_writeRegister(udp_node1, TX_DELAY,  tx_delay);
warplab_writeRegister(udp_node1, TX_LENGTH, tx_length);
warplab_writeRegister(udp_node1, TX_MODE,   tx_mode);
if (channel > 0)
    warplab_setRadioParameter(udp_node1, CARRIER_CHANNEL, channel);
    warplab_setRadioParameter(udp_node2, CARRIER_CHANNEL, channel);
end

% First transmitter: 
if n_tx_ext > 0
    % Radio 1 is the first radio in WARPLab 4x4: 
    warplab_setRadioParameter(udp_node1, RADIO1_TXGAINS, ... 
        (tx1_gain_rf + tx1_gain_bb*2^16));
    warplab_writeSMWO(udp_node1, RADIO1_TXDATA, ...
                      tx_block_pa_filt_up(:,1).');
    warplab_sendCmd(udp_node1, RADIO1_TXEN, packetNum);
    warplab_sendCmd(udp_node1, RADIO1TXBUFF_TXEN, packetNum);
    if print_debug
        fprintf('length(tx_block_pa_filt_up(:,1).'') = %d', ...
                 length(tx_block_pa_filt_up(:,1).'));
    end
end

% Second transmitter: 
if n_tx_ext > 1
    % Radio 2 is the second radio in WARPLab 4x4: 
    warplab_setRadioParameter(udp_node1, RADIO2_TXGAINS, ... 
        (tx2_gain_rf + tx2_gain_bb*2^16));
    warplab_writeSMWO(udp_node1, RADIO2_TXDATA, ...
                      tx_block_pa_filt_up(:,2).');
    warplab_sendCmd(udp_node1, RADIO2_TXEN, packetNum);
    warplab_sendCmd(udp_node1, RADIO2TXBUFF_TXEN, packetNum);
end

% Third transmitter: 
if n_tx_ext > 2
    % Radio 3 is the third radio in WARPLab 4x4: 
    warplab_setRadioParameter(udp_node1, RADIO3_TXGAINS, ... 
        (tx3_gain_rf + tx3_gain_bb*2^16));
    warplab_writeSMWO(udp_node1, RADIO3_TXDATA, ...
                      tx_block_pa_filt_up(:,3).');
    warplab_sendCmd(udp_node1, RADIO3_TXEN, packetNum);
    warplab_sendCmd(udp_node1, RADIO3TXBUFF_TXEN, packetNum);
end

% Fourth transmitter: 
if n_tx_ext > 3
    % Radio 4 is the fourth radio in WARPLab 4x4: 
    warplab_setRadioParameter(udp_node1, RADIO4_TXGAINS, ... 
        (tx4_gain_rf + tx4_gain_bb*2^16));
    warplab_writeSMWO(udp_node1, RADIO4_TXDATA, ...
                      tx_block_pa_filt_up(:,4).');
    warplab_sendCmd(udp_node1, RADIO4_TXEN, packetNum);
    warplab_sendCmd(udp_node1, RADIO4TXBUFF_TXEN, packetNum);
end

% Globally for all transmitters on node 1: 
if n_tx_ext > 0
    warplab_sendCmd(udp_node1, TX_START, packetNum);
end

% First receiver: 
if n_rx > 0
    % Radio 1 is the first radio in WARPLab 4x4: 
    warplab_setRadioParameter(udp_node2, RADIO1_RXGAINS, ...
        (rx1_gain_bb + rx1_gain_rf*2^16));
    warplab_sendCmd(udp_node2, RADIO1_RXEN, packetNum);
    warplab_sendCmd(udp_node2, RADIO1RXBUFF_RXEN, packetNum);
end

% Second receiver: 
if n_rx > 1
    % Radio 2 is the second radio in WARPLab 4x4: 
    warplab_setRadioParameter(udp_node2, RADIO2_RXGAINS, ...
        (rx2_gain_bb + rx2_gain_rf*2^16));
    warplab_sendCmd(udp_node2, RADIO2_RXEN, packetNum);
    warplab_sendCmd(udp_node2, RADIO2RXBUFF_RXEN, packetNum);
end

% Third receiver: 
if n_rx > 2
    % Radio 3 is the third radio in WARPLab 4x4: 
    warplab_setRadioParameter(udp_node2, RADIO3_RXGAINS, ...
        (rx3_gain_bb + rx3_gain_rf*2^16));
    warplab_sendCmd(udp_node2, RADIO3_RXEN, packetNum);
    warplab_sendCmd(udp_node2, RADIO3RXBUFF_RXEN, packetNum);
end

% Fourth receiver: 
if n_rx > 3
    % Radio 4 is the fourth radio in WARPLab 4x4: 
    warplab_setRadioParameter(udp_node2, RADIO4_RXGAINS, ...
        (rx4_gain_bb + rx4_gain_rf*2^16));
    warplab_sendCmd(udp_node2, RADIO4_RXEN, packetNum);
    warplab_sendCmd(udp_node2, RADIO4RXBUFF_RXEN, packetNum);
end

% Globally for all receivers on node 2: 
if n_rx > 0
    warplab_sendCmd(udp_node2, RX_START, packetNum);
end

% Start transmission: 
warplab_sendSync(udp_sync);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the received samples and disable radio resources
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First transmitter: 
if n_tx_ext > 0
    warplab_sendCmd(udp_node1, RADIO1TXBUFF_TXDIS, packetNum);
    warplab_sendCmd(udp_node1, RADIO1_TXDIS, packetNum);
end

% Second transmitter: 
if n_tx_ext > 1
    warplab_sendCmd(udp_node1, RADIO2TXBUFF_TXDIS, packetNum);
    warplab_sendCmd(udp_node1, RADIO2_TXDIS, packetNum);
end

% Third transmitter: 
if n_tx_ext > 2
    warplab_sendCmd(udp_node1, RADIO3TXBUFF_TXDIS, packetNum);
    warplab_sendCmd(udp_node1, RADIO3_TXDIS, packetNum);
end

% Fourth transmitter: 
if n_tx_ext > 3
    warplab_sendCmd(udp_node1, RADIO4TXBUFF_TXDIS, packetNum);
    warplab_sendCmd(udp_node1, RADIO4_TXDIS, packetNum);
end

% First receiver: 
if n_rx > 0
    [rx1_data_raw] = warplab_readSMRO(udp_node2, RADIO1_RXDATA, ... 
                                      tx_length + tx_delay);
    [rx1_data, rx1_RxOTR] = warplab_processRawRxData(rx1_data_raw);
    rx1_data_len = length(rx1_data);
    [rx1_rssi_raw] = warplab_readSMRO(udp_node2, RADIO1_RSSIDATA, ...
                                      ceil((tx_length + tx_delay) / 8));
    [rx1_rssi] = warplab_processRawRSSIData(rx1_rssi_raw);
    warplab_sendCmd(udp_node2, RADIO1RXBUFF_RXDIS, packetNum);
    warplab_sendCmd(udp_node2, RADIO1_RXDIS, packetNum);
end

% Second receiver: 
if n_rx > 1
    [rx2_data_raw] = warplab_readSMRO(udp_node2, RADIO2_RXDATA, ... 
                                      tx_length + tx_delay);
    [rx2_data, rx2_RxOTR] = warplab_processRawRxData(rx2_data_raw);
    rx2_data_len = length(rx2_data);
    [rx2_rssi_raw] = warplab_readSMRO(udp_node2, RADIO2_RSSIDATA, ...
                                      ceil((tx_length + tx_delay) / 8));
    [rx2_rssi] = warplab_processRawRSSIData(rx2_rssi_raw);
    warplab_sendCmd(udp_node2, RADIO2RXBUFF_RXDIS, packetNum);
    warplab_sendCmd(udp_node2, RADIO2_RXDIS, packetNum);
end

% Third receiver: 
if n_rx > 2
    [rx3_data_raw] = warplab_readSMRO(udp_node2, RADIO3_RXDATA, ... 
                                      tx_length + tx_delay);
    [rx3_data, rx3_RxOTR] = warplab_processRawRxData(rx3_data_raw);
    rx3_data_len = length(rx3_data);
    [rx3_rssi_raw] = warplab_readSMRO(udp_node2, RADIO3_RSSIDATA, ...
                                      ceil((tx_length + tx_delay) / 8));
    [rx3_rssi] = warplab_processRawRSSIData(rx3_rssi_raw);
    warplab_sendCmd(udp_node2, RADIO3RXBUFF_RXDIS, packetNum);
    warplab_sendCmd(udp_node2, RADIO3_RXDIS, packetNum);
end

% Fourth receiver: 
if n_rx > 3
    [rx4_data_raw] = warplab_readSMRO(udp_node2, RADIO4_RXDATA, ... 
                                      tx_length + tx_delay);
    [rx4_data, rx4_RxOTR] = warplab_processRawRxData(rx4_data_raw);
    rx4_data_len = length(rx4_data);
    [rx4_rssi_raw] = warplab_readSMRO(udp_node2, RADIO4_RSSIDATA, ...
                                      ceil((tx_length + tx_delay) / 8));
    [rx4_rssi] = warplab_processRawRSSIData(rx4_rssi_raw);
    warplab_sendCmd(udp_node2, RADIO4RXBUFF_RXDIS, packetNum);
    warplab_sendCmd(udp_node2, RADIO4_RXDIS, packetNum);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute and plot the FFT(s) of the received signal(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_spectrum
    figure
    if n_rx > 0
        rx1_fft_len = 2^nextpow2(rx1_data_len);
        rx1_fft = fftshift(fft(rx1_data, rx1_fft_len)/rx1_data_len);
        rx1_fft_f = f_sam/2 * linspace(-1, 1, rx1_fft_len);
        subplot(1, n_rx, 1);
        plot(rx1_fft_f/10^6, abs(rx1_fft));
        title('Spectrum at RX1');
        xlabel('Frequency (MHz)');
        ylabel('Magnitude');
        xlim([-12 12]);
    end
    if n_rx > 1
        rx2_fft_len = 2^nextpow2(rx2_data_len);
        rx2_fft = fftshift(fft(rx2_data, rx2_fft_len)/rx2_data_len);
        rx2_fft_f = f_sam/2 * linspace(-1, 1, rx2_fft_len);
        subplot(1, n_rx, 2);
        plot(rx2_fft_f/10^6, abs(rx2_fft));
        title('Spectrum at RX2');
        xlabel('Frequency (MHz)');
        ylabel('Magnitude');
        xlim([-12 12]);
    end
    if n_rx > 2
        rx3_fft_len = 2^nextpow2(rx3_data_len);
        rx3_fft = fftshift(fft(rx3_data, rx3_fft_len)/rx3_data_len);
        rx3_fft_f = f_sam/2 * linspace(-1, 1, rx3_fft_len);
        subplot(1, n_rx, 3);
        plot(rx3_fft_f/10^6, abs(rx3_fft));
        title('Spectrum at RX3');
        xlabel('Frequency (MHz)');
        ylabel('Magnitude');
        xlim([-12 12]);
    end
    if n_rx > 3
        rx4_fft_len = 2^nextpow2(rx4_data_len);
        rx4_fft = fftshift(fft(rx4_data, rx4_fft_len)/rx4_data_len);
        rx4_fft_f = f_sam/2 * linspace(-1, 1, rx4_fft_len);
        subplot(1, n_rx, 4);
        plot(rx4_fft_f/10^6, abs(rx4_fft));
        title('Spectrum at RX4');
        xlabel('Frequency (MHz)');
        ylabel('Magnitude');
        xlim([-12 12]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the waveform(s) of the received signal(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_waveform
    figure
    if n_rx > 0
        subplot(2, n_rx, 1);
        plot(real(rx1_data));
        title('Receiver 1 I');
        xlabel('Sample'); 
        ylabel('Amplitude');
        axis([0 2^14 -1 1]);
        subplot(2, n_rx, n_rx+1);
        plot(imag(rx1_data));
        title('Receiver 1 Q');
        xlabel('Sample'); 
        ylabel('Amplitude');
        axis([0 2^14 -1 1]);
    end
    if n_rx > 1
        subplot(2, n_rx, 2);
        plot(real(rx2_data));
        title('Receiver 2 I');
        xlabel('Sample'); 
        ylabel('Amplitude');
        axis([0 2^14 -1 1]);
        subplot(2, n_rx, n_rx+2);
        plot(imag(rx2_data));
        title('Receiver 2 Q');
        xlabel('Sample'); 
        ylabel('Amplitude');
        axis([0 2^14 -1 1]);
    end
    if n_rx > 2
        subplot(2, n_rx, 3);
        plot(real(rx3_data));
        title('Receiver 3 I');
        xlabel('Sample'); 
        ylabel('Amplitude');
        axis([0 2^14 -1 1]);
        subplot(2, n_rx, n_rx+3);
        plot(imag(rx3_data));
        title('Receiver 3 Q');
        xlabel('Sample'); 
        ylabel('Amplitude');
        axis([0 2^14 -1 1]);
    end
    if n_rx > 3
        subplot(2, n_rx, 4);
        plot(real(rx4_data));
        title('Receiver 4 I');
        xlabel('Sample'); 
        ylabel('Amplitude');
        axis([0 2^14 -1 1]);
        subplot(2, n_rx, n_rx+4);
        plot(imag(rx4_data));
        title('Receiver 4 Q');
        xlabel('Sample'); 
        ylabel('Amplitude');
        axis([0 2^14 -1 1]);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the Received Signal Strength Indicator (RSSI in dBm) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The conversion is based on the following radio and RSSI vaue
% specifications: 
% rx_gain_rf = 3: rssi_avg=0 is -100dBm; rssi_avg=1023 is -30dBm.
% rx_gain_rf = 2: rssi_avg=0 is  -85dBm; rssi_avg=1023 is -15dBm.
% rx_gain_rf = 1: rssi_avg=0 is  -70dBm; rssi_avg=1023 is   0dBm. 

if n_rx > 0
    rx1_rssi_avg = mean(rx1_rssi);
    rx1_rssi_avg_dBm = (70/1023) * rx1_rssi_avg - 70 - (rx1_gain_rf-1)*15;
    fprintf('\nRX1 RSSI dBm = %5.2f\n',rx1_rssi_avg_dBm);
end
if n_rx > 1
    rx2_rssi_avg = mean(rx2_rssi);
    rx2_rssi_avg_dBm = (70/1023) * rx2_rssi_avg - 70 - (rx2_gain_rf-1)*15;
    fprintf('RX2 RSSI dBm = %5.2f\n',rx2_rssi_avg_dBm);
end
if n_rx > 2
    rx3_rssi_avg = mean(rx3_rssi);
    rx3_rssi_avg_dBm = (70/1023) * rx3_rssi_avg - 70 - (rx3_gain_rf-1)*15;
    fprintf('RX3 RSSI dBm = %5.2f\n',rx3_rssi_avg_dBm);
end
if n_rx > 3
    rx4_rssi_avg = mean(rx4_rssi);
    rx4_rssi_avg_dBm = (70/1023) * rx4_rssi_avg - 70 - (rx4_gain_rf-1)*15;
    fprintf('RX4 RSSI dBm = %5.2f\n',rx4_rssi_avg_dBm);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Downconvert from int_f to baseband
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n_rx > 0
    rx1_time = t_sam * (0 : 1 : rx1_data_len - 1);
    rx1_data_bb = rx1_data .* exp(-j * 2 * pi * int_f * rx1_time);
    rx1_data_bb = rx1_data_bb.';
end
if n_rx > 1
    rx2_time = t_sam * (0 : 1 : rx2_data_len - 1);
    rx2_data_bb = rx2_data .* exp(-j * 2 * pi * int_f * rx2_time);
    rx2_data_bb = rx2_data_bb.';
end
if n_rx > 2
    rx3_time = t_sam * (0 : 1 : rx3_data_len - 1);
    rx3_data_bb = rx3_data .* exp(-j * 2 * pi * int_f * rx3_time);
    rx3_data_bb = rx3_data_bb.';
end
if n_rx > 3
    rx4_time = t_sam * (0 : 1 : rx4_data_len - 1);
    rx4_data_bb = rx4_data .* exp(-j * 2 * pi * int_f * rx4_time);
    rx4_data_bb = rx4_data_bb.';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter the received signal with a matched filter (matched to 
% the SRRC pulse shaping filter) and detect the preamble.
% Then downsample the baseband signal if it was oversampled.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n_rx > 0
    rx1_data_bb_mf = rcosflt(rx1_data_bb, 1, filter_samp, ...
                             'Fs/filter', rrcfilter);
    % Correlate with the reference matrix to find preamble sequence 
    correlation = abs((rx1_data_bb_mf(1:corr_window).')*reference_matrix);
    preamble_start = find(correlation == max(correlation)); 
    first_sample_index = preamble_start + preamble_us_len; 
    % Allow to manually fine tune the first sample position: 
    first_sample_index = first_sample_index + samp_offset;
    if oversamp > 1
        rx1_data_bb = rx1_data_bb_mf(first_sample_index:end);
        rx1_data_bb = downsample(rx1_data_bb, oversamp);
    else
        rx1_data_bb = rx1_data_bb(first_sample_index:end);
    end

    % Truncate the receive vector if too long:
    if length(rx1_data_bb) > size(tx_block, 1)
        rx1_data_bb = rx1_data_bb(1:size(tx_block, 1));
    end
    % Pad the receive vector with zeros if too short:
    if length(rx1_data_bb) < size(tx_block, 1)
        rx1_data_bb(length(rx1_data_bb) + 1 : size(tx_block, 1)) = 0;
    end
    rx_block(:,1) = rx1_data_bb;
end

if n_rx > 1
    if oversamp > 1
        rx2_data_bb_mf = rcosflt(rx2_data_bb, 1, filter_samp, ...
                                 'Fs/filter', rrcfilter);
        rx2_data_bb = rx2_data_bb_mf(first_sample_index:end);
        rx2_data_bb = downsample(rx2_data_bb, oversamp);
    else
        rx2_data_bb = rx2_data_bb(first_sample_index:end);
    end

    % Truncate the receive vector if too long:
    if length(rx2_data_bb) > size(tx_block, 1)
        rx2_data_bb = rx2_data_bb(1:size(tx_block, 1));
    end
    % Pad the receive vector with zeros if too short:
    if length(rx2_data_bb) < size(tx_block, 1)
        rx2_data_bb(length(rx2_data_bb) + 1 : size(tx_block, 1)) = 0;
    end
    rx_block(:,2) = rx2_data_bb;
end

if n_rx > 2
    if oversamp > 1
        rx3_data_bb_mf = rcosflt(rx3_data_bb, 1, filter_samp, ...
                                 'Fs/filter', rrcfilter);
        rx3_data_bb = rx3_data_bb_mf(first_sample_index:end);
        rx3_data_bb = downsample(rx3_data_bb, oversamp);
    else
        rx3_data_bb = rx3_data_bb(first_sample_index:end);
    end

    % Truncate the receive vector if too long:
    if length(rx3_data_bb) > size(tx_block, 1)
        rx3_data_bb = rx3_data_bb(1:size(tx_block, 1));
    end
    % Pad the receive vector with zeros if too short:
    if length(rx3_data_bb) < size(tx_block, 1)
        rx3_data_bb(length(rx3_data_bb) + 1 : size(tx_block, 1)) = 0;
    end
    rx_block(:,3) = rx3_data_bb;
end

if n_rx > 3
    if oversamp > 1
        rx4_data_bb_mf = rcosflt(rx4_data_bb, 1, filter_samp, ...
                                 'Fs/filter', rrcfilter);
        rx4_data_bb = rx4_data_bb_mf(first_sample_index:end);
        rx4_data_bb = downsample(rx4_data_bb, oversamp);
    else
        rx4_data_bb = rx4_data_bb(first_sample_index:end);
    end

    % Truncate the receive vector if too long:
    if length(rx4_data_bb) > size(tx_block, 1)
        rx4_data_bb = rx4_data_bb(1:size(tx_block, 1));
    end
    % Pad the receive vector with zeros if too short:
    if length(rx4_data_bb) < size(tx_block, 1)
        rx4_data_bb(length(rx4_data_bb) + 1 : size(tx_block, 1)) = 0;
    end
    rx_block(:,4) = rx4_data_bb;
end

% Close sockets
pnet('closeall');
