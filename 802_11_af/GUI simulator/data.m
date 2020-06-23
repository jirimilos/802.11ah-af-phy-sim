function [ DATA ] = data(PSDU,SIGNAL,T)
% Data function creates DATA field of the PPDU IEEE 802.11
% PSDU [bits] - input user data of the PPDU
% SIGNAL [bits] - input SIGNAL field of the PPDU containing RATE and LENGTH
% T[table] - input table containing MCS parameters

% check PSDU field
if( length(PSDU) == 0)
    error('PSDU is empty!');
end

% check SIGNAL field
if( length(SIGNAL) != 24)
    error('SIGNAL field is not correct! Have to be 24 bits length!');
end

LENGTH = bi2di(SIGNAL(6:17)); % extract LENGTH of DATA from SIGNAL field

RATE=SIGNAL(1:4); % extract RATE of DATA from SIGNAL field

% check data rate for 20MHz channel spacing
switch(RATE)
    case [1 1 0 1]
        DataRate = '1.1'; % 6/5.625 Mb/s
    case [1 1 1 1]
        DataRate = '1.6'; % 9/5.625 Mb/s
    case [0 1 0 1]
        DataRate = '2.1'; % 12/5.625 Mb/s
    case [0 1 1 1]
        DataRate = '3.2'; % 18/5.625 Mb/s
    case [1 0 0 1]
        DataRate = '4.3'; % 24/5.625 Mb/s
    case [1 0 1 1]
        DataRate = '6.4'; % 36/5.625 Mb/s       
    case [0 0 0 1]
        DataRate = '8.5'; % 48/5.625 Mb/s
    case [0 0 1 1]
        DataRate = '9.6'; % 54/5.625 Mb/s  
    otherwise error('RATE of SIGNAL field incorrect!')
end

% read number of data bits per OFDM symbol for DataRate from SIGNAL field
N_DBPS = T(DataRate,:).N_DBPS;

% define SERVICE bits
scrambler_init = zeros(1,7); % scrambler initialization
service_reserved = zeros(1,9) % reserved bits
SERVICE = [scramler_init service_reserved];

TAIL = zeros(1,6); % define TAIL bits

%define PAD bits
N_SYM = ceil((16+8*LENGTH+6)/N_DBPS); % calculate number of OFDM symbols
N_DATA = N_SYM * N_DBPS; % calculate number of bits in the DATA field
N_PAD = N_DATA-(16+8*LENGTH+6); % calculate number of pad bits
PAD = zeros(1,N_PAD); % define PAD bits of size N_PAD

DATA = [SERVICE PSDU TAIL PAD]; % define output DATA of PPDU

