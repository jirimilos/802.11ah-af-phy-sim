function [ I_filt, Q_filt ] = RRC_filter( ofdm_input, upsampling, beta, span)
% Function filters OFDM signal using Root Raised Cosine filter and returns
% filtered I and Q parts of OFDM signal
% ofdm_input - OFDM signal
% upsampling - upsample coeficient
% beta - RRC filter roll off factor
% span - RRC filter span

% split ofdm signal to real and imag parts
t_Real = real(ofdm_input); 
t_Imag = imag(ofdm_input);

t_Real_up = upsample(t_Real, upsampling);
t_Imag_up = upsample(t_Imag, upsampling);

% RRC filter parameters
sps = upsampling; % Samples per symbol (oversampling factor)

% raised cosine FIR pulse-shaping filter
RRC = rcosdesign(beta,span,sps,'sqrt'); % generate RRC filter coefficients

% set filter delay
t_Real_up(end+span*sps/2)=0;
t_Imag_up(end+span*sps/2)=0;

% apply RRC filter to real and imag parts
I_filt = filter(RRC,1,t_Real_up);
Q_filt = filter(RRC,1,t_Imag_up);

% filter delay compensating
I_filt = I_filt(1,(sps/2)*span+1 : end);
Q_filt = Q_filt(1,(sps/2)*span+1 : end);

end