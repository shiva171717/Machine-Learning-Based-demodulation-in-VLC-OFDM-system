clc;
close all;
clear all;

% Read CSV data
received_bits = csvread('7_44_all.csv'); % Assuming received_bits is 780000x4
Bit_sequence = csvread('Bit_sequence.csv'); % Assuming Bit_sequence is 40000x1
len = length(Bit_sequence) / 4; % Assuming QAM16 log2(16) = 4

% Remove the first row from received_bits
received_bits(1, :) = [];

% Initialize BER matrix
BER = zeros(6, 13); % 6 bandwidths, 13 SNR values

% Loop over each 130000x4 segment to calculate BER
for seg = 1:6
    % Extract the segment corresponding to each bandwidth (130000x4)
    segment_start = (seg-1) * 130000 + 1;
    segment_end = seg * 130000;
    recvd_bits_segment_full = received_bits(segment_start:segment_end, :);

    % Loop over each 10000x4 segment within the 130000x4 segment
    for i = 1:13
        % Extract the i-th 10000x4 segment from received_bits
        sub_segment_start = (i-1) * len + 1;
        sub_segment_end = i * len;
        recvd_bits_segment = recvd_bits_segment_full(sub_segment_start:sub_segment_end, :);

        % Reshape the segment from 10000x4 to 1x40000
        recvd_bits_segment = reshape(recvd_bits_segment', 1, []);

        % Reshape Bit_sequence from 40000x1 to 1x40000
        Bit_sequence_reshaped = reshape(Bit_sequence, 1, []);

        % Calculate the number of errors for this segment
        num_errors = sum(Bit_sequence_reshaped ~= recvd_bits_segment);

        % Calculate BER for this segment, ensuring no zeros
        BER(seg, i) = num_errors / length(Bit_sequence_reshaped);
        if BER(seg, i) == 0
            BER(seg, i) = eps; % Add a small value to avoid zero BER
        end
    end
end

% Define SNR values
SNR = 0:5:60;

% Plot SNR vs BER using semilogy for each bandwidth
figure;
hold on;
bandwidths = [6, 7, 9, 11, 13, 15]; % MHz

for seg = 1:6
    semilogy(SNR, BER(seg, :), '-o', 'DisplayName', sprintf('%d MHz', bandwidths(seg)));
end

xlabel('SNR (dB)');
ylabel('BER');
title('SNR vs BER for Each 130000x4 Bits Segment (Different Bandwidths)');
set(gca, 'YScale', 'log'); % Ensure log scale on y-axis
ylim([1e-6, 1e-1]); % Adjust the y-axis range for better visibility
grid on;
legend('show');
hold off;

display(BER);
