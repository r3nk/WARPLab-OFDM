%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OFDM Framework
%
% Copyright (c) 2012 Robin Klose 
%
% This function initiates an OFDM transmission from a single transmitter
% to a single receiver. It is designed to be used in conjunction with the 
% WARPLab Framework to enable physical transmissions. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OFDM parameter setup
% Provide a structure with the following parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Nrx: [1,4]
%   Number of receivers.
%
% Nsubcar:
%   Number of subcarriers.
%
% Npadl:
%   Number of zero padding bins at low frequencies (should be odd).
%
% Npadh:
%   Number of zero padding bins at high frequencies (should be odd).
%
% Nguard:
%   Number of guard space taps.
%
% Npil: 
%   Number of pilot OFDM symbols in the block header. 
%
% Ndat: 
%   Number of payload/data OFDM symbols in the block.
%
% bps: 
%   Bits per symbol. A corresponding QAM modulation scheme will
%   automatically be used. E.g., for bps = 4, 16-QAM is used. 
%
% cfo_correction: [0,1] 
%   Perform the carrier frequency offset (CFO) correction. The CFO
%   estimation is performed on a sequence of pilot symbols. Thus, a longer
%   sequence of pilots yields a more precise estimation and correction than
%   a short sequence. A blind estimation on data symbols is not performed. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARPLab parameter setup 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% transceiver: 
%   Specifies whether the data is transmitted via WARPLab 2x2, WARPLab 4x4
%   or if the transmission is simulated. Specify a string: 
%   'clean':        Simulation without noise. 
%   'awgn':         Simulation with additive white gaussian noise. 
%   'warplab_2x2':  WARPLab 2x2
%   'warplab_4x4':  WARPLab 4x4
%
% oversamp: 
%   Apply oversampling on data before transmission with WARPLab.
%   1: Disable oversampling
%   Any number > 1: oversampling factor
%
% int_f: Intermediate frequency
%   If Npadl is set to 0, int_f should be > Nc * delta_f / 2.
%   If Npadl is 1 or greater, no subcarriers are placed at 0 Hz and int_f
%   should be set to 0.
%
% carrier_channel: [1,14]
%   Carrier channel used by WARPLab. 
%   Set to 0 to not change a previously setup channel.
%
% samp_offset: [default: 0]
% Manually shift the position of the first sample in the received data
% block. 
% 
% WARPLab TX and RX gain parameters. 
% In warplab_2x2_multi_transceive(), the radios are assigned as follows:
% TX1 is node 1 radio 2. 
% RX1 is node 2 radio 2. RX2 is node 2 radio 3. 
% In warplab_4x4_multi_transceive(), the radios are assigned as follows:
% TX1 is node 1 radio 1. 
% RX1 is node 2 radio 1. RX2 is node 2 radio 2. 
% RX3 is node 2 radio 3. RX4 is node 2 radio 4. 
% 
% tx1_gain_bb: [0, 3] WARPLab TX1 baseband gain.
% tx1_gain_rf: [0,63] WARPLab TX1 radio frequency gain.
% rx1_gain_bb: [0,31] WARPLab RX1 baseband gain.
% rx2_gain_bb: [0,31] WARPLab RX2 baseband gain.
% rx3_gain_bb: [0,31] WARPLab RX3 baseband gain.
% rx4_gain_bb: [0,31] WARPLab RX4 baseband gain.
% rx1_gain_rf: [1, 3] WARPLab RX1 radio frequency gains.
% rx2_gain_rf: [1, 3] WARPLab RX2 radio frequency gains.
% rx3_gain_rf: [1, 3] WARPLab RX3 radio frequency gains.
% rx4_gain_rf: [1, 3] WARPLab RX4 radio frequency gains.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print and plot switches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% print_setup: 
%   Print parameter setup.
%
% print_debug: 
%   Print debug messages. 
%
% print_SNR: 
%   Print measured SNR values. 
%
% plot_filter: 
%   Plot impulse response of pulse shaping filter. 
%   The pulse shaping filter is applied to transmission data only if
%   oversampling is enabled. It is always applied to the synchronization
%   preamble. 
%
% plot_spectrum: 
%   Plot spectrum of received signal.
%
% plot_waveform: 
%   Plot waveform of received signal.
%
% plot_bode: 
%   Plot bode diagram of channel estimate.
%
% plot_SNR: 
%   Plot SNR estimates. 
%
% plot_scatter: 
%   Plot scatter diagram. 
%
% plot_error_lines: 
%   Plot lines for erroneous symbols' origins in scatter plot.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Indices for filenames to save records of transmissions for 
% later evaulation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% name_pref:
%   Filename prefix.
%
% run_major: 
%   Major index of record. 
%   Use this variable to distinguish between records with different
%   parameter settings. 
%
% run_minor:
%   Minor index of record. 
%   Use this variable to distinguish between multiple records with 
%   identical parameter settings. 
% 
% In order to not save the transmission, set both indices to 0.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ofdm_func(params)


% --------------------------------
% Constant system parameters
% --------------------------------

% WARPLab sampling frequency is 40 MHz: 
f_sam = 40e6; % Hz
t_sam = 1 / f_sam; % s


% --------------------------------
% Independent OFDM parameters
% --------------------------------

Nrx     = params.Nrx;
Nsubcar = params.Nsubcar;
Npadl   = params.Npadl;
Npadh   = params.Npadh;
Nguard  = params.Nguard;
Npil    = params.Npil;
Ndat    = params.Ndat;
bps     = params.bps;
cfo_correction = params.cfo_correction;

% --------------------------------
% WARPLab transmission parameters
% --------------------------------

transceiver = params.transceiver; 
oversamp = params.oversamp; 
int_f = params.int_f;
carrier_channel = params.carrier_channel;
samp_offset = params.samp_offset;

% WARPLab TX baseband gain: [0, 3]
tx1_gain_bb = params.tx1_gain_bb;

% WARPLab TX radio frequency gain: [0,63]
tx1_gain_rf = params.tx1_gain_rf;

% WARPLab RX baseband gains: [0,31]
rx1_gain_bb = params.rx1_gain_bb;
rx2_gain_bb = params.rx2_gain_bb;
rx3_gain_bb = params.rx3_gain_bb;
rx4_gain_bb = params.rx4_gain_bb;

% WARPLab RX radio frequency gains: [1, 3]
rx1_gain_rf = params.rx1_gain_rf;
rx2_gain_rf = params.rx2_gain_rf;
rx3_gain_rf = params.rx3_gain_rf;
rx4_gain_rf = params.rx4_gain_rf;

% --------------------------------
% Plot switches
% --------------------------------

print_setup = params.print_setup;
print_debug = params.print_debug;
print_SNR = params.print_SNR;
plot_filter = params.plot_filter;
plot_spectrum = params.plot_spectrum;
plot_waveform = params.plot_waveform;
plot_bode = params.plot_bode;
plot_SNR = params.plot_SNR;
plot_scatter = params.plot_scatter;
plot_error_lines = params.plot_error_lines;


% --------------------------------
% Indices to save measurements:
% --------------------------------

name_pref = params.name_pref;
run_major = params.run_major;
run_minor = params.run_minor;


% --------------------------------
% Channel simulation parameters
% --------------------------------

% Signal to noise ratio in dB for AWGN channel:
awgn_snr_dB = 26; % dB

% Maximum speed of nodes:
v_max = 10; % m/s

% Note: Speed of light is 299792458 m/s
% Speed of wave propagation: 
c = 299792458; % m/s

% Center frequency: 
f_center = 2.4e9; % Hz

% Maximum channel delay (approximately):
tau_max = 400e-9;
% 400ns ~ -20dB in the LTE Extended Pedestrian A model (EPA)


% --------------------------------
% Dependent OFDM parameters
% --------------------------------

% Total number of OFDM subcarriers:
Nc = Nsubcar;

% Total number of OFDM frequency components:
Nt = Nc + Npadl + Npadh;

% Block size: Total number of subsequent OFDM symbols in the block:
Npildat = Npil + Ndat;

% OFDM symbol duration:
t_sym = t_sam * oversamp * Nt; % s
t_sym_cp = t_sam * oversamp * (Nt + Nguard); % s

% OFDM subcarrier spacing:
delta_f = 1 / t_sym; % Hz

% OFDM passband bandwidth:
B_pb = Nc * delta_f; % Hz

% Modulation scheme: see also: help dmod
% modem.qammod, modem.pskmod

% Modulator:
% enmod = modem.qammod(2^bps);
enmod = modem.qammod('M', 2^bps, ...
                     'PhaseOffset', 0, ...
                     'SymbolOrder', 'gray', ...
                     'InputType', 'integer');

% Demodulator: 
demod = modem.qamdemod(enmod);

% --------------------------------
% Dependent variables of simulated channel
% --------------------------------

% Coherence time calculation:
f_doppler_max = (v_max * f_center) / c;
coh_t = 1 / (2*f_doppler_max);
coh_sym = coh_t / t_sym; 

% Coherence bandwidth calculation: 
coh_bw = 1 / tau_max;
coh_scar = coh_bw / delta_f;

% --------------------------------
% Print parameter setup: 
% --------------------------------

if print_setup

fprintf('\n\n===========================================\n');
fprintf('OFDM Parameters\n');
fprintf('===========================================\n');
fprintf('Total # of subcarriers:    %d\n', Nc);
fprintf('Total # of freq. comp.:    %d\n', Nt);
fprintf('Symbol duration:           %.3f µs\n', t_sym * 1e6);
fprintf('Guard spacing:             %d*%d samples\n', oversamp, Nguard);
fprintf('                           %.1f ns\n', ... 
                                    Nguard * t_sam * oversamp * 1e9);
fprintf('Subcarrier spacing:        %.3f kHz\n', delta_f * 1e-3);
fprintf('Lowest frequency:          %.3f kHz\n', ...
        ((Npadl/2) * delta_f) * 1e-3);
fprintf('                           [should be > 30 kHz]\n');
fprintf('Highest frequency:         %.3f MHz\n', ...
        (((Npadl + Nc)/2)) * delta_f * 1e-6);
fprintf('                           [should be < 9.5 MHz]\n');
fprintf('Occupied RF bandwidth:     %.3f MHz\n', B_pb * 1e-6);

if (strcmp(transceiver, 'clean')) || (strcmp(transceiver, 'awgn'))
fprintf('\n\n===========================================\n');
fprintf('Simulated Channel Parameters\n');
fprintf('===========================================\n');
fprintf('Supposed max. mobility:    v_max = %.1f m/s\n', v_max);
fprintf('Coherence time:            coh_t = %.1f µs\n', ...
        coh_t*1e6);
fprintf('                           coh_sym = %.1f symbols\n', ...
        coh_sym);
fprintf('Supposed max. delay:       tau_max = %.1f ns\n', tau_max*1e9);
fprintf('Coherence bandwidth:       coh_bw = %.1f MHz\n', ...
        coh_bw*1e-6);
fprintf('                           coh_scar = %.1f subcarriers\n', ...
        coh_scar);
fprintf('\n');
end
end

% --------------------------------
% Transmitter part
% --------------------------------

% Creation of OFDM transmission vector: A block of OFDM symbols contains
% a sequence of pilot symbols at the beginning, followed by random data. 
% Data generation for pilots and data is not distinguished in the
% transmitter model. The receiver, however, handels symbols either as 
% pilots or as data depending on the symbols' positions in the sequence
% according to the Npil and Ndat parameters. 

% Matrix organization: 
% The first dimensions codes the data within an OFDM symbol. 
% The second dimension holds sequential OFDM symbols (pilots or tones).

% Generate random values. Subsequently, these values are either used
% as pilot symbols or as data depending on their position in the matrix.
tx_val = randint(Nc, Npildat, 2^bps);
% Nc#Npildat

% Map values to symbols in constellation plane: 
tx_sym = enmod.modulate(tx_val);
% Nc#Npildat

% Maximum I/Q value of symbols: 
tx_sym_max_i  = max(abs(real(tx_sym(:)))); % 1#1
tx_sym_max_q  = max(abs(imag(tx_sym(:)))); % 1#1
tx_sym_max_iq = max(tx_sym_max_i, tx_sym_max_q); % 1#1


% Allocate a zero padded matrix. This matrix also contains frequency
% components that are not used as subcarriers by the transmitter.
% Frequency bins that are used as subcarriers are filled with values
% subsequently.
tx_sym_full_pad = zeros(Nt, Npildat);
% Nt#Npildat

% Positions to fill the padded matrix with values: 
pad_low_idx = 1 + ceil(Npadl / 2);
pad_low_len = floor(Nc / 2);
pad_hig_idx = pad_low_idx + pad_low_len + Npadh;
pad_hig_len = ceil(Nc / 2);

% Fill padded matrix with values: 
tx_sym_full_pad(pad_low_idx : pad_low_idx + pad_low_len - 1, :) = ...
                tx_sym(1 : pad_low_len, :);
tx_sym_full_pad(pad_hig_idx : pad_hig_idx + pad_hig_len - 1, :) = ...
                tx_sym(pad_low_len + 1 : end, :);

% Perform IFFT on the data to obtain OFDM symbols.
% The IFFT is calculated along columns by default. 
% The iDFT is per se normalized by 1/Nt, so multiply with sqrt(Nt)
% in order to obtain the OFDM normalization of 1/sqrt(Nt): 
tx_sym_ofdm = ifft(tx_sym_full_pad(:,:)) * sqrt(Nt);
% Nt#Npildat

% Add cyclic prefix of length Nguard: 
if Nguard > 0
    tx_sym_ofdm_cp = [tx_sym_ofdm(end-Nguard+1:end,:); tx_sym_ofdm];
else
    tx_sym_ofdm_cp = tx_sym_ofdm;
end
% (Nt+Nguard)#Npildat

% Transform matrix to column vector containing all OFDM symbols: 
tx_ofdm_vector = reshape(tx_sym_ofdm_cp, [], 1);
% (Nt+Nguard)*Npildat#1

% Maximum I/Q value of OFDM vector:
tx_ofdm_vector_max_iq = max([abs(real(tx_ofdm_vector)); ...
                             abs(imag(tx_ofdm_vector))]); 
% 1#1

% Normalize I/Q samples to [-1;1] range to meet WARP's DAC requirements:
tx_ofdm_vector = tx_ofdm_vector / tx_ofdm_vector_max_iq;


% --------------------------------
% Transmission
% --------------------------------

% Transmission without noise: 
if (strcmp(transceiver, 'clean'))
    rx_ofdm_vectors = repmat(tx_ofdm_vector, [1 Nrx]);
end

% Transmission with AWGN: 
if (strcmp(transceiver, 'awgn'))
    rx_ofdm_vector_noisy = awgn(tx_ofdm_vector, awgn_snr_dB);
    rx_ofdm_vectors = repmat(rx_ofdm_vector_noisy, [1 Nrx]);
end

% Set up WARPLab transceive parameters: 
transceive.n_tx = 1;
transceive.n_rx = Nrx;
transceive.tx_block = tx_ofdm_vector;
transceive.oversamp = oversamp;
transceive.int_f = int_f;
transceive.channel = carrier_channel;
transceive.samp_offset = samp_offset;
transceive.tx1_gain_bb = tx1_gain_bb;
transceive.tx2_gain_bb = 0;
transceive.tx3_gain_bb = 0;
transceive.tx4_gain_bb = 0;
transceive.tx1_gain_rf = tx1_gain_rf;
transceive.tx2_gain_rf = 0;
transceive.tx3_gain_rf = 0;
transceive.tx4_gain_rf = 0;
transceive.rx1_gain_bb = rx1_gain_bb;
transceive.rx2_gain_bb = rx2_gain_bb;
transceive.rx3_gain_bb = rx3_gain_bb;
transceive.rx4_gain_bb = rx4_gain_bb;
transceive.rx1_gain_rf = rx1_gain_rf;
transceive.rx2_gain_rf = rx2_gain_rf;
transceive.rx3_gain_rf = rx3_gain_rf;
transceive.rx4_gain_rf = rx4_gain_rf;
transceive.print_debug = print_debug;
transceive.plot_filter = plot_filter;
transceive.plot_spectrum = plot_spectrum;
transceive.plot_waveform = plot_waveform;

% Transmission via WARPLab 2x2: 
if (strcmp(transceiver, 'warplab_2x2'))
    rx_ofdm_vectors = warplab_2x2_multi_transceive(transceive);
end

% Transmission via WARPLab 4x4: 
if (strcmp(transceiver, 'warplab_4x4'))
    rx_ofdm_vectors = warplab_4x4_multi_transceive(transceive);
end

fprintf('TX dim: %dx%d\t RX dim: %dx%d\n', ...
        size(tx_ofdm_vector, 1), size(tx_ofdm_vector, 2), ...
        size(rx_ofdm_vectors, 1), size(rx_ofdm_vectors, 2));

% --------------------------------
% Receiver part
% --------------------------------

% Truncation is handled internally by warplab_nxn_transmit, so the 
% receive vector starts with the first sample of the OFDM transmission. 

% Each column is an OFDM symbol with a cyclic prefix:
rx_sym_ofdm_cp = reshape(rx_ofdm_vectors, Nt+Nguard, Npildat, Nrx);

% Remove cyclic prefix:
rx_sym_ofdm = rx_sym_ofdm_cp(1 + Nguard : end, :, :);
% Nt#Npildat#Nrx

% Allocate matrices for OFDM symbols: 
rx_sym_full_pad = zeros(Nt, Npildat, Nrx);

% Perform FFT on the received data to retain the complex symbols.
% The FFT is calculated along columns by default. Devide by
% sqrt(Nt) to retain the OFDM normalization of 1/sqrt(Nt): 
for rx = 1:Nrx
    rx_sym_full_pad(:,:,rx) = fft(rx_sym_ofdm(:,:,rx)) / sqrt(Nt);
end

% Extract I/Q symbols from zero padded matrices: 
rx_sym_raw = [rx_sym_full_pad(pad_low_idx : ... 
                              pad_low_idx + pad_low_len - 1, :, :); ...
              rx_sym_full_pad(pad_hig_idx : ...
                              pad_hig_idx + pad_hig_len - 1, :, :)];
% Nc#Npildat#Nrx


% ----------------
% CFO estimation and correction:
% ----------------

if (cfo_correction)
    % This CFO estimation algorithm averages the phase drift over a
    % sequence of successive pilot symbols. The phase drift between
    % successive symbols in the sequence must not exceed +/- pi. 

    % Each column contains an estimate of the phase shift of its
    % corresponding symbol. However, the received symbols are subjected
    % to noise as well which has an impact on CFO estimation. 
    % Noise has the same impact on every symbol while the CFO affects
    % symbols accumulatively, i.e., the n-th OFDM symbol has a phase shift
    % of (n-1) times the phase offset between successive symbols.
    cfo_sym_phase_shifts = zeros(Nc, Npil, Nrx);
    for rx = 1:Nrx
        cfo_sym_phase_shifts(:, :, rx) = ... 
            angle(rx_sym_raw(:, 1:Npil, rx) ./ tx_sym(:, 1:Npil));
    end
    % in range of ]-pi, pi]

    % Calculate the differences of phase shifts of successive symbols:
    cfo_sym_phase_diffs = diff(cfo_sym_phase_shifts, 1, 2);
    % Nc#(Npil-1)#Nrx in range of ]-2pi, 2pi[
    
    % The sum of the square of the phase shifts is used as a measure for
    % the precision of the CFO estimate on a particular subcarrier. Some
    % subcarriers may be exposed to channel impairments which might
    % impact the CFO estimation. The higher the phase differences between
    % symbols, the higher is the impact on CFO estimation. 
    cfo_sym_phase_diffs_2sum = sum(cfo_sym_phase_diffs.^2, 2);
    cfo_sym_phase_diffs_weights = ...
                        max(cfo_sym_phase_diffs_2sum(:)) ./ ... 
                            cfo_sym_phase_diffs_2sum;
    % Nc#1#Nrx
    
    % Averaging the phase differences between successive symbols reduces
    % the impact of noise and exposes the systematic error given by the
    % continuous phase drift per symbol. The arithmetic mean, however, 
    % is not effective for angles. Therefore, take the argument of the 
    % mean of the complex numbers with the angles as argument: 
    cfo_sym_drift_mean_rx_subcar = ... 
                      angle(sum(exp(j*cfo_sym_phase_diffs), 2));
    % Nc#1#Nrx

    % The CFO affects all subcarriers equally. However, some subcarriers
    % may be subjected to channel impairments and fading und thus induce 
    % errors on the CFO estimate. Therefore, average the CFO drift over all
    % subcarriers after weighting with their channel quality measure:
    cfo_sym_drift_mean_rx = ... 
                      angle(sum(cfo_sym_phase_diffs_weights .* ...
                          exp(j*cfo_sym_drift_mean_rx_subcar), 1));
    % 1#1#Nrx

    % The CFO affects all receivers equally because they are driven by the
    % same clock on the same WARP board. Therefore, average the CFO drift
    % estimates over all receivers:
    cfo_sym_drift_mean = angle(sum(exp(j*cfo_sym_drift_mean_rx)));

    % The frequency difference is the phase drift per symbol 
    % divided by the symbol duration including the guard space: 
    cfo_delta_f = cfo_sym_drift_mean / (t_sym_cp * 2 * pi); % Hz

    fprintf('delta f = %.1f Hz\n', cfo_delta_f);

    % Estimated accumulated phase drift for pilot and data symbols: 
    cfo_sym_drift_est_acc = (0:Npildat-1) * cfo_sym_drift_mean;

    % Compensate accumulated phase drift: 
    rx_sym_raw = bsxfun(@times, rx_sym_raw, ...
                        exp(-j * cfo_sym_drift_est_acc));

end


% ----------------
% Estimate channel:
% ----------------

% Estimate channel coefficients for each subcarrier and each receiver: 
H_est = zeros(Nc, Npil, Nrx);
for rx = 1:Nrx
    H_est(:, :, rx) = rx_sym_raw(:, 1:Npil, rx) ./ tx_sym(:, 1:Npil);
end

% Mean of all pilots per subcarrier:
% H_est_mean_subcar contains a channel estimate for each individual
% subcarrier. 
H_est_mean_subcar = sum(H_est, 2) / Npil;
% Nc#1#Nrx

% Maximum absolute channel gain for subcarrier estimate:
H_est_mean_subcar_abs_max = max(abs(H_est_mean_subcar)); 
% 1#1#Nrx


% ----------------
% Pilot-aided SNR estimation:
% ----------------

% The noise power corresponds to the variance of the noise.
% Interpolate the square of the absolute distance of the measured 
% constellation points to their corresponding constellation point over 
% a sequence of pilot symbols to estimate the noise variance: 

% Extract pilot symbols from the received data blocks:
rx_sym_raw_pil = rx_sym_raw(:, 1:Npil, :);
% Nc#Npil#Nrx

% Estimate the origin of pilots with channel estimates for subcarriers: 
rx_pil_est_subcar = bsxfun(@rdivide, rx_sym_raw_pil, H_est_mean_subcar);
% Nc#Npil#Nrx

% Square of absolute distance between received pilots and their origin: 
noise_est_abs_dist_square = zeros(Nc, Npil, Nrx);
for rx = 1:Nrx
    noise_est_abs_dist_square(:,:,rx) = ... 
        abs(rx_pil_est_subcar(:,:,rx) - tx_sym(:, 1:Npil)).^2;
end
% Nc#Npil#Nrx

% Estimation of simga square (noise power): 
noise_est_sigma_square = sum(noise_est_abs_dist_square, 2) / Npil;
% Nc#1#Nrx

% Estimate the SNR: 
SNR_est_pil_subcar = ((tx_sym_max_i)^2 + (tx_sym_max_q)^2) ./ ...
                        noise_est_sigma_square; 
SNR_est_pil_subcar_dB = 10 * log10(SNR_est_pil_subcar);
% Nc#1#Nrx

SNR_est_pil_subcar_dB_max_rx = max(SNR_est_pil_subcar_dB, [], 1);
%1#1#Nrx

SNR_est_pil_subcar_dB_max = max(SNR_est_pil_subcar_dB_max_rx(:));
% 1#1#1

% Mean SNR values of receivers: 
SNR_est_pil_rx_mean = mean(SNR_est_pil_subcar, 1);
SNR_est_pil_rx_mean_dB = 10 * log10(SNR_est_pil_rx_mean);
% 1#1#Nrx

if print_SNR
    for rx = 1:Nrx
        fprintf('Mean estimated SNR RX%d: %3.1f dB\n', ... 
                             rx, SNR_est_pil_rx_mean_dB(rx));
    end
end

% ----------------
% Apply channel estimate:
% ----------------

% Remove pilot symbols from the received block:
rx_sym_raw_dat = rx_sym_raw(:, 1 + Npil : end, :);
% Nc#Ndat#Nrx

% Estimate the received symbols with channel estimates for subcarriers: 
rx_sym_dat = bsxfun(@rdivide, rx_sym_raw_dat, H_est_mean_subcar);
% Nc#Ndat#Nrx

% Maximum I/Q value of received symbols per receiver: 
rx_sym_dat_vec = reshape(rx_sym_dat, [], Nrx); % Nc*Ndat#Nrx
rx_sym_dat_rx_max_i = max(abs(real(rx_sym_dat_vec)), [], 1); % 1#Nrx
rx_sym_dat_rx_max_q = max(abs(imag(rx_sym_dat_vec)), [], 1); % 1#Nrx
rx_sym_dat_rx_max_iq = max([rx_sym_dat_rx_max_i; ... 
                            rx_sym_dat_rx_max_q], [], 1); % 1#Nrx
% rx_sym_dat_max_iq = max(rx_sym_dat_rx_max_iq(:)); % 1


% ----------------
% Demodulate symbols:
% ----------------

rx_val_dat = zeros(Nc, Ndat, Nrx);
for rx = 1:Nrx
    rx_val_dat(:,:,rx) = demod.demodulate(rx_sym_dat(:,:,rx));
end


% ----------------
% Error statistics:
% ----------------

% Get linear error indices for individual receivers:
sym_err_rx(Nrx).lin = [];

% Symbol and bit errors at receivers:
sym_err_rx_cnt = zeros(Nrx);
sym_err_rx_rat = zeros(Nrx);
bit_err_rx_cnt = zeros(Nrx);
bit_err_rx_rat = zeros(Nrx);

for rx = 1:Nrx

    % Linear error indices for individual receivers: 
    sym_err_rx(rx).lin = ...
        find(tx_val(:, 1+Npil:end) - rx_val_dat(:,:,rx) ~= 0);

    % Number and rate of symbol errors: 
    [sec, ser] = symerr(tx_val(:, 1+Npil:end), ... 
                        rx_val_dat(:,:,rx));
    sym_err_rx_cnt(rx) = sec;
    sym_err_rx_rat(rx) = ser;

    % Number and rate of bit errors: 
    [bec, ber] = biterr(tx_val(:, 1+Npil:end), ... 
                        rx_val_dat(:,:,rx), bps);
    bit_err_rx_cnt(rx) = bec;
    bit_err_rx_rat(rx) = ber;
end

% Symbol and bit errors for individual subcarriers:
sym_err_scar_cnt = zeros(Nc, Nrx);
sym_err_scar_rat = zeros(Nc, Nrx);
bit_err_scar_cnt = zeros(Nc, Nrx);
bit_err_scar_rat = zeros(Nc, Nrx);
for rx = 1:Nrx
    for scar = 1:Nc
        [sec, ser] = symerr(tx_val(scar, 1+Npil:end), ...
                            rx_val_dat(scar, :, rx));
        sym_err_scar_cnt(scar, rx) = sec;
        sym_err_scar_rat(scar, rx) = ser;
        
        [bec, ber] = biterr(tx_val(scar, 1+Npil:end), ... 
                            rx_val_dat(scar, :, rx), bps);
        bit_err_scar_cnt(scar, rx) = bec;
        bit_err_scar_rat(scar, rx) = ber;
    end
end

% Print error statistics:
for rx = 1:Nrx
    fprintf('RX%d:   ', rx);
    fprintf('SEC:%d   SER:%02.3f%%   BEC:%d   BER:%02.3f%%\n', ... 
            sym_err_rx_cnt(rx), 100*sym_err_rx_rat(rx), ... 
            bit_err_rx_cnt(rx), 100*bit_err_rx_rat(rx));
end


% --------------------------------
% Plots
% --------------------------------

% Plot Bode diagram of measured channel over subcarriers: 
if plot_bode
    figure;
    set(gcf,'Color',[1, 1, 1]);
    for rx = 1:Nrx
        subplot(2, Nrx, rx);
        plot(abs(H_est_mean_subcar(:,:,rx)), ':.');
        axis([1, Nc, 0, 1.1*H_est_mean_subcar_abs_max(rx)]);
        xlabel('Subcarrier');
        ylabel('Magnitude');
        title(sprintf('H(f) at RX%d', rx));
        subplot(2, Nrx, Nrx + rx);
        plot(angle(H_est_mean_subcar(:,:,rx)), ':.');
        axis([1, Nsubcar, -pi, pi]);
        set(gca,'YTick',-pi:pi/2:pi);
        set(gca,'YTickLabel',{'-pi','-pi/2','0','pi/2','pi'});
        xlabel('Subcarrier');
        ylabel('Phase');
    end
    % Show only one legend: 
    legend(cellstr('TX'));
end


% SNR over subchannels: 
if plot_SNR
    figure;
    set(gcf,'Color',[1, 1, 1]);
    for rx = 1:Nrx
        subplot(1, Nrx, rx);
        stem(SNR_est_pil_subcar_dB(:,:,rx), ':.');
        if (SNR_est_pil_subcar_dB_max > 0)
            axis([0, Nc+1, 0, 1.1*SNR_est_pil_subcar_dB_max]);
        end
        xlabel('Subcarrier');
        ylabel('dB');
        title(sprintf('SNR at RX%d', rx));
    end
    % Show only one legend: 
    legend(cellstr('TX'));
end


% Scatter plot of received/transmitted constellation points and errors: 
if plot_scatter
    legendstr(Nrx).str = [];
    for rx = 1:Nrx
        % This line gets a new handle for a scatterplot which is afterwards
        % used in the loop. 
        % It is done this way since there is a bug in the Communications
        % Toolbox 3.2 (R14SP3) which makes a scatterplot not appear in a 
        % subplot. See solution ID: 1-23329R 
        h = scatterplot(enmod.modulate(0:1:2^bps-1), 1, 0, 'k*');
        hold on;
        legendstr(rx).str = cellstr('Transmitted');
        
        % Plot received symbols:
        vec_rx_all = reshape(rx_sym_dat(:,:,rx), [], 1);
        scatterplot(vec_rx_all, 1, 0, 'g.', h);
            legendstr(rx).str(2) = ... 
                cellstr(sprintf('TX ~ RX%d ok', rx));
            
        % Plot symbol errors: 
        if (sym_err_rx_cnt(rx) > 0)
            vec_rx = rx_sym_dat(:,:,rx);
            sym_err_rx_symbols = vec_rx(sym_err_rx(rx).lin);
            scatterplot(sym_err_rx_symbols(:), 1, 0, 'r.', h);
            legendstr(rx).str(length(legendstr(rx).str)+1) = ... 
                        cellstr(sprintf('TX ~ RX%d error', rx));
        end

        % Lines to symbol errors:
        if ((sym_err_rx_cnt(rx) > 0) && plot_error_lines)
            vec_tx = tx_sym(:, 1+Npil:end);
            vec_rx = rx_sym_dat(:,:,rx);
            sym_err_tx_symbols = vec_tx(sym_err_rx(rx).lin);
            sym_err_rx_symbols = vec_rx(sym_err_rx(rx).lin);
            for ii = 1:length(sym_err_rx(rx).lin)
                line([real(sym_err_tx_symbols(ii)); ...
                      real(sym_err_rx_symbols(ii))], ...
                     [imag(sym_err_tx_symbols(ii)); ...
                      imag(sym_err_rx_symbols(ii))], ...
                      'linewidth', 0.25, 'color', [1,0,0]);
            end
        end

        % Plot original symbols again in order to be on top:
        scatterplot(enmod.modulate(0:1:2^bps-1), 1, 0, 'k*', h);

        % Show legend:
        legend(legendstr(rx).str);

        title(sprintf('RX%d Constellations', rx));
        max_iq = max(tx_sym_max_iq, rx_sym_dat_rx_max_iq(rx));
        axis([-1.1*max_iq 1.1*max_iq -1.1*max_iq 1.1*max_iq]);
        hold off;
        grid on;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save measurements: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ((run_major > 0) || (run_minor > 0))
    save(sprintf('%s_%04d_%04d', name_pref, run_major, run_minor), ...
            'params', ...
            'f_sam', ... 
            't_sam', ... 
            'awgn_snr_dB', ... 
            'v_max', ... 
            'c', ...
            'f_center', ... 
            'tau_max', ...
            'tx_val', ... 
            'rx_ofdm_vectors', ...
            'cfo_delta_f');
end
