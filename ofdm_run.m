%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OFDM Framework
%
% Copyright (c) 2012 Robin Klose 
%
% This script allows to configure the PHY and visualization parameters 
% and to initiate OFDM transmissions. Transmissions can be recorded to
% files in order to facilitate automated execution in a loop. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

% --------------------------------
% OFDM parameter setup
% --------------------------------

% Number of receiving nodes:
params.Nrx = 4;

% Number of subcarriers:
params.Nsubcar = 64;

% Number of zero padding bins at low frequencies (should be odd):
params.Npadl = 1;

% Number of zero padding bins at high frequencies (should be odd):
params.Npadh = 15;

% Number of guard space taps:
params.Nguard = 8;

% Number of pilot OFDM symbols in the block header:
params.Npil = 20;

% Number of payload/data OFDM symbols in the block:
params.Ndat = 20;

% Bits per symbol:
params.bps = 6;

% Perform a carrier frequency offset (CFO) correction. The CFO estimation
% is performed on a sequence of pilot symbols. Thus, a longer sequence of
% pilots yields a more precise estimation and correction than a short
% sequence. A blind estimation on data symbols is not performed.
params.cfo_correction = 1;

% --------------------------------
% WARPLab transmission parameters
% --------------------------------

% Specifies whether the data is transmitted via WARPLab 2x2, WARPLab 4x4
% or if the transmission is simulated. Specify a string:
% 'clean':        Simulation without noise.
% 'awgn':         Simulation with additive white gaussian noise.
% 'warplab_2x2':  WARPLab 2x2
% 'warplab_4x4':  WARPLab 4x4
params.transceiver = 'warplab_4x4';

% Apply oversampling on data before transmission with WARPLab:
% 1: Disable oversampling
% Any number > 1: oversampling factor
params.oversamp = 2;

% Intermediate frequency:
% If Npadl is set to 0, int_f should be > Nc * delta_f / 2.
% If Npadl is 1 on greater, no subcarriers are placed at 0 Hz and int_f
% should be set to 0.
params.int_f = 0;

% Carrier channel used by WARPLab: [1,14]
% Set to 0 to not change a previously setup channel.
params.carrier_channel = 9;

% Manually fine tune the alignment of the block of received samples.
% Negative values will make the transceiver read in the samples from an
% earlier position in time, positive values cause an artificial delay.
params.samp_offset = 0;

% WARPLab TX and RX gain parameters.
% In warplab_2x2_multi_transceive(), the radios are assigned as follows:
% TX1 is node 1 radio 2. TX2 is node 1 radio 3.
% RX1 is node 2 radio 2. RX2 is node 2 radio 3.
% In warplab_4x4_multi_transceive(), the radios are assigned as follows:
% TX1 is node 1 radio 1. TX2 is node 1 radio 2.
% TX3 is node 1 radio 3. TX4 is node 1 radio 4.
% RX1 is node 2 radio 1. RX2 is node 2 radio 2.
% RX3 is node 2 radio 3. RX4 is node 2 radio 4.

% WARPLab TX baseband gains: [0, 3]
params.tx1_gain_bb = 2;

% WARPLab TX radio frequency gains: [0,63]
params.tx1_gain_rf = 40;

% WARPLab RX baseband gains: [0,31]
params.rx1_gain_bb = 16;
params.rx2_gain_bb = 16;
params.rx3_gain_bb = 16;
params.rx4_gain_bb = 16;

% WARPLab RX radio frequency gains: [1, 3]
params.rx1_gain_rf = 3;
params.rx2_gain_rf = 3;
params.rx3_gain_rf = 3;
params.rx4_gain_rf = 3;


% --------------------------------
% Print and plot switches
% --------------------------------

% Print parameter setup:
params.print_setup = 1;

% Print debug messages:
params.print_debug = 0;

% Print measured SNR values:
params.print_SNR = 1;

% Plot impulse response of pulse shaping filter.
% The pulse shaping filter is applied to transmission data only if
% oversampling is enabled. It is always applied to the synchronization
% preamble.
params.plot_filter = 0;

% Plot spectrum(s) of received signal(s):
params.plot_spectrum = 1;

% Plot waveform(s) of received signal(s):
params.plot_waveform = 1;

% Plot bode diagram of channel estimate:
params.plot_bode = 1;

% Plot SINR estimates:
params.plot_SNR = 1;

% Plot scatter diagram:
params.plot_scatter = 1;

% Plot lines for erroneous symbols' origins in scatter plot:
params.plot_error_lines = 1;

% --------------------------------
% Indices for filenames to save records for later evaluation
% --------------------------------

% Prefix for filename to save experiment:
params.name_pref = 'records/experiment';

% Set both params.run_major and params.run_minor to 0 in order to not save
% the transmission. 

% Major index to identify recorded transmissions: 
params.run_major = 0;

% Minor index to identify recorded transmissions: 
params.run_minor = 0;

% --------------------------------
% Run the transmission(s)
% --------------------------------


ofdm_func(params);

% for ii = 1:32
%     params.run_minor = ii;
%     ofdm_func(params);
% end
