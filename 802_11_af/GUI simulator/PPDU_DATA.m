function [ DATA ] = PPDU_DATA(PSDU,SIGNAL,T)
% Data function creates DATA field of the PPDU IEEE 802.11af
% PSDU [bits] - input user data of the PPDU
% SIGNAL [bits] - input SIGNAL field of the PPDU containing RATE and LENGTH
% T[table] - input table containing MCS parameters

% check PSDU field
if( isempty(PSDU) == 1)
    error('PSDU is empty!');
end

% check SIGNAL field
if( length(SIGNAL) ~= 24)
    error('SIGNAL field is not correct! Have to be 24 bits length!');
end

LENGTH = bi2de(SIGNAL(6:17)); % extract LENGTH of DATA from SIGNAL field

RATE=mat2str(SIGNAL(1:4)); % extract RATE of DATA from SIGNAL field

% check RATE vector in table
T.Properties.RowNames = T.RATE'; % define RATE bits rows

% read number of data bits per OFDM symbol for DataRate from SIGNAL field
N_DBPS = T(RATE,:).N_DBPS; %read coding rate from T

%MCS = MCS{1}; % converts cell to string

% read data rate for RATE from SIGNAL field
%DataRate = T(RATE,:).DataRate;
%DataRate = DataRate{1}; % converts cell to string

%T.Properties.RowNames = T.DataRate'; % define data rate rows
% read number of data bits per OFDM symbol for DataRate from SIGNAL field
%N_DBPS = T(DataRate,:).N_DBPS;

% define SERVICE bits
scrambler_init = zeros(1,7); % scrambler initialization
service_reserved = zeros(1,9); % reserved bits
SERVICE = [scrambler_init service_reserved];

TAIL = zeros(1,6); % define TAIL bits

%define PAD bits
N_SYM = ceil((16+8*LENGTH+6)/N_DBPS); % calculate number of OFDM symbols
N_DATA = N_SYM * N_DBPS; % calculate number of bits in the DATA field
N_PAD = N_DATA-(16+8*LENGTH+6); % calculate number of pad bits
PAD = zeros(1,N_PAD); % define PAD bits of size N_PAD

DATA = [SERVICE PSDU TAIL PAD]; % define output DATA of PPDU
