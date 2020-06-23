function [ DATA ] = PPDU_DATA_TVHT(PSDU, MCS, T)

% check PSDU field
if( isempty(PSDU) == 1)
    error('PSDU is empty!');
end

APEP_LENGTH = length(PSDU)/8; % number of <1;1048575> bytes
N_SERVICE = 16; % number of bits in the SERVICE field
N_TAIL = 6; % number of tail bitrs per BCC encoder
N_ES = 1; % number of BCC endcoders

% check MCS vector in table
T.Properties.RowNames = T.MCS'; % define MCS index rows

% read N_DBPS for selected MCS index
N_DBPS = T(MCS,:).N_DBPS; %read coding rate from T

m_STBC = 1; % 1-> STBC is not used ; 2-> when STBC is used

% calculate number of data symbols in the DATA field
N_SYM = m_STBC * ceil( ( 8 * APEP_LENGTH + N_SERVICE + N_TAIL * N_ES) / (m_STBC * N_DBPS ) );

% calculate PSDU length
PSDU_LENGTH = floor( ( N_SYM * N_DBPS - N_SERVICE - N_TAIL * N_ES ) / 8 );

% define SERVICE bits
scrambler_init = zeros(1,7); % scrambler initialization
service_reserved = 0; % reserved bit
service_CRC = zeros(1,8); % VHT-SIG-B CRC 8 bits

SERVICE = [scrambler_init service_reserved service_CRC];

TAIL = zeros(1,N_TAIL); % define TAIL bits

%define PAD bits
N_PAD = N_SYM * N_DBPS - 8 * PSDU_LENGTH - N_SERVICE - N_TAIL * N_ES;
PAD = zeros(1,N_PAD);

if ( APEP_LENGTH < PSDU_LENGTH )
    N_add_zeros = PSDU_LENGTH - APEP_LENGTH;
    bytes_to_fill_PSDU = zeros(1,N_add_zeros*8);
    PSDU = [PSDU bytes_to_fill_PSDU];
end

DATA = [SERVICE PSDU TAIL PAD]; % define output DATA of PPDU

end

