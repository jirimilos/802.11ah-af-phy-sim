function [ SIGNAL ] = signal( MCS, data_length)
% Signal function creates SIGNAL field of the PPDU IEEE 802.11af
% data_rate [Mb/s] - input data rate for DATA field of the PPDU
% data_length [B] - input data length for DATA field of the PPDU
% SIGNAL [bits] - output SIGNAL field of the PPDU IEEE 802.11af

% check selected MCS index
% X = 5.625 for 8 MHz channel | X = 7.5 for 6 or 7 MHz channel
switch(MCS)
    case '0'
        RATE = [1 1 0 1]; % 6/X Mb/s
    case '1'
        RATE = [1 1 1 1]; % 9/X Mb/s
    case '2'
        RATE = [0 1 0 1]; % 12/X Mb/s
    case '3'
        RATE = [0 1 1 1]; % 18/X Mb/s
    case '4'
        RATE = [1 0 0 1]; % 24/X Mb/s
    case '5'
        RATE = [1 0 1 1]; % 36/X Mb/s       
    case '6'
        RATE = [0 0 0 1]; % 48/X Mb/s
    case '7'
        RATE = [0 0 1 1]; % 54/X Mb/s  
    otherwise error('MCS index incorrect!')
end

% check data length
if (data_length > 4095)
    error('Data length incorrect! Possible data length in range <0;4095> bytes')
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