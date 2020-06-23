function [ ofdm_output, ofdm_output_with_pilots ] = ofdm_demodulator2( ofdm_input, FORMAT, CH_BANDWIDTH, BCU_WIDTH, cyclic_prefix, extra_zeros_flag, equalizer_flag)
%   Function process OFDM demodulation (removes GI interval) for OFDM symbols
%   and returns complex subcarriers output
%   ofdm_input - [complex numbers] data input to OFDM demodulatorB
%   FORMAT - [string] PDDU format used 
%   CH_BANDWIDTH - [string] channel bandwidth used
%   BCU_WIDTH - [string] basic channel unit width
%   ofdm_output - [complex numbers] data output from OFDM demodulator
%   cyclic_prefix - [string] flag for short/normal CP
%   extra_zeros_flag - [int] flag for extra zeros added in ofdm modulator
%   (1/0)
%   equalizer_flag - [int] flag if equalizer used (1/0)

if nargin<7
  equalizer_flag = -1;
end

switch(FORMAT)
    case 'non-HT' % if non-HT format choosen NON_HT duplicate modulation used
        N_SD = 96; % data subcarriers
        N_SP = 8; % pilot subcarriers
        N_ST = N_SD + N_SP; % 104 subcarriers in total per BCU
        
    case 'VHT' % if VHT format choosen
        N_SD = 108; % data subcarriers
        N_SP = 6; % pilot subcarriers
        N_ST = N_SD + N_SP; % 114 subcarriers in total per BCU
    otherwise error('Incorrect format used! Could be non-HT/VHT')
end

switch(BCU_WIDTH)
    case {'6 MHz','8 MHz'}
        BCU_W = 144; % number of subcarriers per one BCU
    case '7 MHz' 
        BCU_W = 168; % number of subcarriers per one BCU
    otherwise error('Incorrect basic channel unit width used! Could be 6 MHz/7 MHz/8 MHz')
end

switch(CH_BANDWIDTH)
    case 'TVHT_W'      
        IFFT_size = 128;
        switch(FORMAT)
            case 'non-HT'
                subcarrier_non_HT_indexes = [-26:-1,1:26]; % subcarrier indexes per non-HT mode
                subcarriers_DC = zeros(1,(10+1)); % subcarriers (1 DC + 10 NULLs)

                subcarriers_non_HT = [ % subcarriers for non-HT mode
                    subcarrier_non_HT_indexes(find(subcarrier_non_HT_indexes == -26):find(subcarrier_non_HT_indexes == -1)) 0 ... 
                    subcarrier_non_HT_indexes(find(subcarrier_non_HT_indexes == 1):find(subcarrier_non_HT_indexes == 26))];

                subcarriers_non_HT_dup = [ % subcarriers for non-HT duplicate mode
                    subcarriers_non_HT subcarriers_DC subcarriers_non_HT];

                % NULLs addition to fill IFFT block 
                fft_nulls =IFFT_size - length(subcarriers_non_HT_dup); % number of NULLs to fill IFFT block
                a_nulls = ceil(fft_nulls/2); % NULLs on the begin
                b_nulls = fft_nulls-a_nulls; % NULLs on the end
                subcarriers_non_HT_dup = [zeros(1,a_nulls) subcarriers_non_HT_dup zeros(1,b_nulls)]; % IFFT block mapping

                pilot1_indexes = find( subcarriers_non_HT_dup == -21 );
                pilot2_indexes = find( subcarriers_non_HT_dup == -7 );
                pilot3_indexes = find( subcarriers_non_HT_dup == 7 );
                pilot4_indexes = find( subcarriers_non_HT_dup == 21 );

                pilot_indexes = [pilot1_indexes ... % subcarriers with index -21
                                pilot2_indexes ... % subcarriers with index -7
                                pilot3_indexes ... % subcarriers with index +7
                                pilot4_indexes];% subcarriers with index +21

                % DC nulls mapping
                nulls_indexes = find( subcarriers_non_HT_dup == 0 );
            
            case 'VHT'

                subcarrier_VHT_indexes = [-58:-2,2:58]; % subcarrier indexes per VHT mode
                subcarriers_DC = zeros(1,(2+1)); % subcarriers (1 DC + 2 NULLs)

                subcarriers_VHT = [ % subcarriers for VHT mode include DC
                    subcarrier_VHT_indexes(find(subcarrier_VHT_indexes == -58):find(subcarrier_VHT_indexes == -2)) subcarriers_DC ... 
                    subcarrier_VHT_indexes(find(subcarrier_VHT_indexes == 2):find(subcarrier_VHT_indexes == 58))];

                % NULLs addition to fill IFFT block 
                fft_nulls =IFFT_size - length(subcarriers_VHT); % number of NULLs to fill IFFT block
                a_nulls = ceil(fft_nulls/2); % NULLs on the begin
                b_nulls = fft_nulls-a_nulls; % NULLs on the end
                subcarriers_VHT = [zeros(1,a_nulls) subcarriers_VHT zeros(1,b_nulls)]; % IFFT block mapping

                % pilots mapping

                pilot1_indexes = find( subcarriers_VHT == -53 );
                pilot2_indexes = find( subcarriers_VHT == -25 );
                pilot3_indexes = find( subcarriers_VHT == -11 );
                pilot4_indexes = find( subcarriers_VHT == 11 );
                pilot5_indexes = find( subcarriers_VHT == 25 );
                pilot6_indexes = find( subcarriers_VHT == 53 );

                pilot_indexes = [pilot1_indexes ... % subcarriers with index -53
                                pilot2_indexes ... % subcarriers with index -25
                                pilot3_indexes ... % subcarriers with index -11
                                pilot4_indexes ...% subcarriers with index +11
                                pilot5_indexes ... % subcarriers with index +25
                                pilot6_indexes]; % subcarriers with index +53

                % DC nulls mapping
                nulls_indexes = find( subcarriers_VHT == 0 );
        end        
        
    case 'TVHT_2W' 
        IFFT_size = 512;
    case 'TVHT_4W'
     IFFT_size = 1024;
    case 'TVHT_W+W'
        IFFT_size = 128;
    case 'TVHT_2W+2W' 
        IFFT_size = 512;
    otherwise error('Incorrect channel bandwidth used! Could be TVHT_(W/2W/4W/W+W/2W+2W)')
end

switch(cyclic_prefix)
    case 'GI'
        T_GI = floor(BCU_W/4); % GI duration
    case 'GIS' 
        T_GI = floor(BCU_W/8); % GI short duration
    otherwise error('Incorrect cyclic prefix! Could be GI/GIS')
end


%% S2P - Serial to Parallel conversion

N_S = BCU_W + T_GI;
N_SYM = length(ofdm_input)/N_S;
S2P = zeros(N_SYM,N_S); 

k=0;

for i=1:N_SYM
    for j=1:N_S
        S2P(i,j) = ofdm_input(k*N_S+j);
    end
    k = k + 1;
end

%% Remove CP

% removing cyclic prefix
CP_start = T_GI + 1; % start position for OFDM symbol
subcarriers_per_BCU = S2P(:,CP_start:end); % OFDM symbols without cyclic prefix

% extract IFFT part from BCU (removing effective guard band)
guard_band_nulls = BCU_W - IFFT_size;
subcarriers = subcarriers_per_BCU(:,((guard_band_nulls/2)+1):(end-(guard_band_nulls/2)));

% apply FFT operation on subcarriers
for n = 1:N_SYM
    subcarriers(n,:) = fft(subcarriers(n,:),IFFT_size);
end

%% Subcarriers demapping

subcarriers = fftshift(subcarriers,2); % shifted IFFT block mapping

% data demapping

indexes = [pilot_indexes nulls_indexes]; %(not used due to equalizer)

%subcarriers(:,nulls_indexes) = []; % remove all null subcarriers

if ( equalizer_flag == 1 )
    subcarriers = equalizerZF(subcarriers, pilot_indexes, CH_BANDWIDTH, FORMAT);
end

%% P2S conversion for output include pilots and nulls

ofdm_output_with_pilots = [];

subcarriers2 = ifftshift(subcarriers,2); % shifted IFFT block mapping

% apply FFT operation on subcarriers
for n = 1:N_SYM
    subcarriers2(n,:) = ifft(subcarriers2(n,:),IFFT_size);
end

if( strcmp(FORMAT,'non-HT') )
    k=0;
    for i=1:N_SYM % convert parallel stream to serial
        for j=1:(IFFT_size)
            if ( i == N_SYM && j > ((IFFT_size)/2) && extra_zeros_flag == 1)
                break;
            else
                ofdm_output_with_pilots = [ofdm_output_with_pilots subcarriers2(i,j)];            
            end
        end
        k = k + 1;
    end
elseif( strcmp(FORMAT,'VHT') )
    for i=1:N_SYM % convert parallel stream to serial
        for j=1:IFFT_size
            ofdm_output_with_pilots = [ofdm_output_with_pilots subcarriers2(i,j)];            
        end
    end
end

subcarriers(:,indexes) = []; % remove all pilot+zeros subcarriers

a=mean( abs( subcarriers( end, 1 : ( N_SD / 2 ) ) ) );
b=mean( abs( subcarriers( end, ( ( N_SD / 2 ) + 1 ) : end) ) );

% if ( a > b  && abs(a-b) > 0.2)
%     extra_zeros_detected = 1; % flag is set to extra zeros detected
% else
%     extra_zeros_detected = 0; % flag is not set to extra zeros detected
% end

% if ( sum( fix( abs( subcarriers(end,((N_SD/2)+1):end) ) ) ) == 0)
%     extra_zeros_detected = 1; % flag is set to extra zeros detected
% else
%     extra_zeros_detected = 0; % flag is not set to extra zeros detected
% end


%% P2S conversion

ofdm_output = [];

if( strcmp(FORMAT,'non-HT') )
    k=0;
    for i=1:N_SYM % convert parallel stream to serial
        for j=1:N_SD
            if ( i == N_SYM && j > (N_SD/2) && extra_zeros_flag == 1)
                break;
            else
                ofdm_output = [ofdm_output subcarriers(i,j)];            
            end
        end
        k = k + 1;
    end
elseif( strcmp(FORMAT,'VHT') )
    for i=1:N_SYM % convert parallel stream to serial
        for j=1:N_SD
            ofdm_output = [ofdm_output subcarriers(i,j)];            
        end
    end
end
