function [ SIGNAL ] = signal( data_rate, data_length)
% Signal function creates SIGNAL field of the PPDU IEEE 802.11
% data_rate [Mb/s] - input data rate for DATA field of the PPDU
% data_length [B] - input data length for DATA field of the PPDU
% SIGNAL [bits] - output SIGNAL field of the PPDU IEEE 802.11

% check selected data rate for 20MHz channel spacing
switch(data_rate)
    case '6'
        RATE = [1 1 0 1]; % 6 Mb/s
    case '9'
        RATE = [1 1 1 1]; % 9 Mb/s
    case '12'
        RATE = [0 1 0 1]; % 12 Mb/s
    case '18'
        RATE = [0 1 1 1]; % 18 Mb/s
    case '24'
        RATE = [1 0 0 1]; % 24 Mb/s
    case '36'
        RATE = [1 0 1 1]; % 36 Mb/s       
    case '48'
        RATE = [0 0 0 1]; % 48 Mb/s
    case '54'
        RATE = [0 0 1 1]; % 54 Mb/s  
    otherwise error('Data rate incorrect! Possible data rates 6/9/12/18/24/36/48/54 Mb/s')
end

% check data length
if (data_length > 65535)
    error('Data length incorrect! Possible data length in range <0;65535> bytes')
else
    LENGTH = de2bi(data_length,12); % transform to 12-bit LENGTH field
end

%define reserved bit of SIGNAL field
R = 0;

%calculate parity bit (even parity)
bits = [RATE R LENGTH]; % define 0 - 16 bits to calculate parity bit
P=mod(sum(bits),2); % define parity bit

%define tail bits
TAIL = zeros(1,6); % define 6 zero tail bits

SIGNAL = [bits P TAIL]; % define SIGNAL field of PPDU frame 