clc;
close all;
clear all;

% Parameters
active_subc = 500;
data_block = 40;
M = 16;
Eb = 1;
block_size = 1024;
cp_len = 8;
EbNo_dB = 44;
No = Eb * 10.^(-1 * EbNo_dB / 10);


%% you can change data rate and cutoff frequency of LPF 
data_rate   =   15e6; 
Tb  =   1/data_rate;
filter_order = 1;
pass_band = 7e6;  % cutoff freq of LPF 
n_samp=10;  
samp_pos = 1;           % sample position from center of symbol


% Generate Bit Sequence
Bit_sequence = pnseq11((active_subc * 2) * data_block);
csvwrite('Bit_sequence.csv', Bit_sequence.');
% Map bits to symbols
data_source = zeros(1, length(Bit_sequence) / log2(M));
for i = 1:log2(M):length(Bit_sequence)
    bit_str = Bit_sequence(i:i+log2(M)-1);
    data_source((i-1)/log2(M) + 1) = bin2dec(num2str(bit_str));
end

% QAM Modulation
qQAM_modulated_data = qammod(data_source, M);

% Display QAM constellation
scatterplot(qQAM_modulated_data);
title('QAM Modulated Transmitted Data');

% OFDM Modulation
num_cols = length(qQAM_modulated_data) / active_subc;
data_matrix = reshape(qQAM_modulated_data, active_subc, num_cols);

% Add zero padding and Hermitian symmetry
ifft_data_matrix = zeros(block_size, num_cols);
ifft_data_matrix(2:active_subc+1, :) = data_matrix;            
ifft_data_matrix(end-active_subc+1:end, :) = flipud(conj(data_matrix)); 

% Perform IFFT 
ifft_data = ifft(ifft_data_matrix, block_size,1);
ifft_data_cp = [ifft_data(end-cp_len+1:end, :); ifft_data];

% Serialize OFDM signal
ofdm_signal = ifft_data_cp(:).';


% Plot OFDM Signal
figure;
plot(real(ofdm_signal));
xlabel('Time');
ylabel('Amplitude');
title('OFDM Signal');
grid on;

no_of_sample = length(ofdm_signal);     % number of samples
Total_time_window = no_of_sample/data_rate;
delta_T = Total_time_window/no_of_sample; 
sampling_freq = n_samp/delta_T;

over_sampled_OFDM = repelem(ofdm_signal,n_samp);
% BER Simulation
BER = 0;

    snr = EbNo_dB;
    No = 1/(10^(snr/10));
    num_errors = 0;
    total_bits = 0;

    filtered=lowpass_bessel_filter_output(filter_order,2*pi*pass_band,over_sampled_OFDM,sampling_freq);

    % Add AWGN noise
    noise = sqrt(No/2) * (randn(size(filtered)) + 1i * randn(size(filtered)));
    received_signal = filtered + noise;
    downsamp_sig=received_signal((n_samp/2)+samp_pos:n_samp:end);  % sampling at end of bit period


    % Reshape to matrix form
    recvd_signal_matrix = reshape(downsamp_sig, block_size + cp_len, num_cols);

    % Remove cyclic prefix
    recvd_signal_matrix(1:cp_len, :) = [];

    % Perform FFT
    fft_data_matrix = fft(recvd_signal_matrix, block_size, 1);
    comp_data = fft_data_matrix(2:active_subc+1, :); % Extract data subcarriers

    % Convert to serial stream
    recvd_serial_data = comp_data(:).';
    recvd_serial_data = recvd_serial_data;
    scatterplot(recvd_serial_data);

    % QAM Demodulation
    qQAM_demodulated_data = qamdemod(recvd_serial_data, M);

    % Map symbols back to bits
    received_bits = de2bi(qQAM_demodulated_data, log2(M), 'left-msb')';
    received_bits = received_bits(:)';

    % Calculate BER
    num_errors = num_errors + sum(Bit_sequence ~= received_bits);
    total_bits = total_bits + length(Bit_sequence);

    BER = num_errors / total_bits;
    display(BER);




csvwrite('7_44recvd_serial_data.csv', recvd_serial_data.');
csvwrite('7_44data_source.csv', data_source.');
