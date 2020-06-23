function [ output_signal, carrier_I, carrier_Q ] = IQ_modulator(I, Q, upsampling, fc, fs)
% Function process IQ modulation for I and Q parts of filtered OFDM signal
% I - In-phase part of filtered OFDM signal
% Q - Quadrature part of filtered OFDM signal
% fs - sampling frequency
% fc - frequency carrier

% Upconverter

%fs = (8e+6)/144; % one symbol frequency X kHz

%upsampling = 10;

delta_f = fs * upsampling; % X MHz

t = 0:(1/delta_f):(length(I)-1)/delta_f;

% create real and imag carriers
carrier_I = cos(2*pi*fc*t);
carrier_Q = sin(2*pi*fc*t);

% IQ modulation
signal_I = I.*carrier_I;
signal_Q = Q.*carrier_Q;

output_signal = signal_I + signal_Q;

end

