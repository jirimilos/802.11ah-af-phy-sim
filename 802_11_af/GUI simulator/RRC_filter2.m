function [ output_signal ] = RRC_filter2( I, Q, downsampling, beta, span)
% Function filters IQ demodulated signal using Root Raised Cosine filter and returns
% OFDM signal
% I - Real part of IQ demodulated signal
% Q - Imag part of IQ demodulated signal
% downsapling - downsampling frequency (should be equal to usampling frequency)
% beta - RRC filter roll off factor
% span - RRC filter span

% RRC filter parameters
sps = downsampling; % Samples per symbol (oversampling factor)

% raised cosine FIR pulse-shaping filter
RRC = rcosdesign(beta,span,sps,'sqrt'); % generate RRC filter coefficients

% set filter delay
I(end+span*sps/2)=0;
Q(end+span*sps/2)=0;

% filter I and Q signals
I_filt = filter(RRC,0.5,I);
Q_filt = filter(RRC,0.5,Q);

% filter delay compensating
I_compensated = I_filt(1,(sps/2)*span+1 : end);
Q_compensated = Q_filt(1,(sps/2)*span+1 : end);

% downsample I and Q filtered signals

I_received = downsample(I_compensated, downsampling);
Q_received = downsample(Q_compensated, downsampling);

% sum I and Q received signals
output_signal = I_received + 1i*Q_received;

end

