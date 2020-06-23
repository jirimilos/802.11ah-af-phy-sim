function [ ofdm_output, extra_zeros_flag ] = ofdm_modulator2( ofdm_input, FORMAT, CH_BANDWIDTH, BCU_WIDTH, cyclic_prefix)
%   Function process OFDM modulation of complex modulated subcarriers
%   and products OFDM symbols to output

%   ofdm_input - [complex numbers] data input to OFDM modulator
%   FORMAT - [string] PDDU format used 
%   CH_BANDWIDTH - [string] channel bandwidth used
%   BCU_WIDTH - [string] basic channel unit width
%   ofdm_output - [complex numbers] data output from OFDM modulator

switch(FORMAT)
    case 'non-HT' % if non-HT format choosen NON_HT duplicate modulation used
        
        N_SD = 96; % data subcarriers
        N_SP = 8; % pilot subcarriers
        N_ST = N_SD + N_SP; % 104 subcarriers in total per BCU
        
        subcarrier_non_HT_indexes = [-26:-1,1:26]; % subcarrier indexes per non-HT mode
        subcarriers_DC = zeros(1,(10+1)); % subcarriers (1 DC + 10 NULLs)

        subcarriers_non_HT = [ % subcarriers for non-HT mode
            subcarrier_non_HT_indexes(find(subcarrier_non_HT_indexes == -26):find(subcarrier_non_HT_indexes == -1)) 0 ... 
            subcarrier_non_HT_indexes(find(subcarrier_non_HT_indexes == 1):find(subcarrier_non_HT_indexes == 26))];

        subcarriers_non_HT_dup = [ % subcarriers for non-HT duplicate mode
            subcarriers_non_HT subcarriers_DC subcarriers_non_HT];
        
    case 'VHT' % if VHT format choosen
        
        N_SD = 108; % data subcarriers
        N_SP = 6; % pilot subcarriers
        N_ST = N_SD + N_SP; % 114 subcarriers in total per BCU
        
        subcarrier_VHT_indexes = [-58:-2,2:58]; % subcarrier indexes per VHT mode
        subcarriers_DC = zeros(1,(2+1)); % subcarriers (1 DC + 2 NULLs)

        subcarriers_VHT = [ % subcarriers for VHT mode include DC
            subcarrier_VHT_indexes(find(subcarrier_VHT_indexes == -58):find(subcarrier_VHT_indexes == -2)) subcarriers_DC ... 
            subcarrier_VHT_indexes(find(subcarrier_VHT_indexes == 2):find(subcarrier_VHT_indexes == 58))];
                
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
        N_SEG = 1; % number of channel segments
    case 'TVHT_2W' 
        IFFT_size = 512;
        N_SEG = 2; % number of channel segments
        switch(FORMAT)
            case 'non-HT'
            case 'VHT'
                subcarriers_NULLs = zeros(1,(13*2+1)); % subcarriers (X NULLs)
                subcarriers_VHT = [subcarriers_VHT subcarriers_NULLs subcarriers_VHT];
        end
    case 'TVHT_4W'
        IFFT_size = 1024;
        N_SEG = 4; % number of channel segments
    case 'TVHT_W+W'
        IFFT_size = 128;
        N_SEG = 2; % number of channel segments
        switch(FORMAT)
            case 'non-HT'
            case 'VHT'
                subcarriers_NULLs = zeros(1,(13*2+1)); % subcarriers (X NULLs)
                subcarriers_VHT = [subcarriers_VHT subcarriers_NULLs subcarriers_VHT];
        end
    case 'TVHT_2W+2W' 
        IFFT_size = 512;
        N_SEG = 4; % number of channel segments
    otherwise error('Incorrect channel bandwidth used! Could be TVHT_(W/2W/4W/W+W/2W+2W)')
end

switch(FORMAT)
    case 'non-HT'

        % NULLs addition to fill IFFT block 
        fft_nulls =IFFT_size - length(subcarriers_non_HT_dup); % number of NULLs to fill IFFT block
        a_nulls = ceil(fft_nulls/2); % NULLs on the begin
        b_nulls = fft_nulls-a_nulls; % NULLs on the end
        subcarriers_non_HT_dup = [zeros(1,a_nulls) subcarriers_non_HT_dup zeros(1,b_nulls)]; % IFFT block mapping

        % pilots mapping

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

switch(cyclic_prefix)
    case 'GI'
        T_GI = floor(BCU_W/4); % GI duration
    case 'GIS' 
        T_GI = floor(BCU_W/8); % GI short duration
    otherwise error('Incorrect cyclic prefix! Could be GI/GIS')
end

N_SD = N_SD * N_SEG;
N_SYM = length(ofdm_input)/N_SD; % number of OFDM symbols per BCU
    
if( strcmp(FORMAT,'non-HT') )

    % calculations for non-HT duplicate mode
    N_SYM_NON_HT = length(ofdm_input)/(N_SD/2); % number of OFDM symbols in non-HT mode

    if ( mod( N_SYM_NON_HT , 2 ) == 1 ) % if odd number of OFDM symbols
        extra_zeros_flag = 1; % flag is set 1 because of extra nulls
        ofdm_input = [ofdm_input zeros(1,(N_SD/2))]; % correction by adding 48 nulls as SIGNAL/DATA part
        N_SYM = ceil(N_SYM); % number of fft blocks correction
    else
        extra_zeros_flag = 0; % flag is set 0 because of not extra nulls
    end

    %% pilots initialization
    % pn sequence of pilots polarity for 127 OFDM symbols
    p127 = [1,1,1,1, -1,-1,-1,1, -1,-1,-1,-1, 1,1,-1,1, -1,-1,1,1, -1,1,1,-1, 1,1,1,1, 1,1,-1,1, 1,1,-1,1, 1,-1,-1,1, 1,1,-1,1, -1,-1,-1,1, -1,1,-1,-1, 1,-1,-1,1, 1,1,1,1, -1,-1,1,1, -1,-1,1,-1, 1,-1,1,1, -1,-1,-1,1, 1,-1,-1,-1, -1,1,-1,-1, 1,-1,1,1, 1,1,-1,1, -1,1,-1,1, -1,-1,-1,-1, -1,1,-1,1, 1,-1,1,-1, 1,1,1,-1, -1,1,-1,-1, -1,1,1,1, -1,-1,-1,-1, -1,-1,-1];
    % pilots sequence of pilots -21, -7, +7, +21
    P = [1 1 1 -1];
    p127_coeff = ceil(N_SYM/127);
    if ( p127_coeff == 1)
        pn = p127;
    else % extend pn sequency according to number of OFDM symbols
        pn = repmat(p127,1,p127_coeff);
    end
    
elseif( strcmp(FORMAT,'VHT') )
    
    %% pilots initialization
    % pilots sequence of pilots -53 -25 -11 +11 +25 +53 for 1 spatial stream
    P = [1, 1, 1, -1, -1, 1];
    
    % pilots values for N_SYM symbols  
    pilots = zeros(N_SYM, N_SP);
    
    for i = 1:N_SYM
        
        %calculate P sequence indexes for corresponding N_SYM symbol
        psi_1 = mod(i, N_SP);
        psi_2 = mod(i+1, N_SP);
        psi_3 = mod(i+2, N_SP);
        psi_4 = mod(i+3, N_SP);
        psi_5 = mod(i+4, N_SP);
        psi_6 = mod(i+5, N_SP);
        
        %pn sequence for corresponding N_SYM symbol
        pn = [P(psi_1+1) P(psi_2+1) P(psi_3+1) P(psi_4+1) P(psi_5+1) P(psi_6+1)];
        
        %fill pilots with pn sequence
        pilots(i,:) = pn;
        
    end
    
end
    
% % pilots sequence of pilots -53, -39, -26, -12, +12, +26, +39, +53
% PSI = [1 1 1 -1 -1 1 1 1]; % 8 pilots vector
% 
% pilots = zeros(N_SYM,8); % all pilots N_SYM x 8
% 
% % calculate pilots for all OFDM symbols
% for i=1:N_SYM
%     for j=1:8
%         index = mod( ( ( i - 1 ) + ( j - 1 ) ) , 8 ) + 1;
%         pilots(i,j) = PSI(index);
%     end
% end

% 
% if flag == 0 % SIGNAL part case
%     
%     
% 
% else % DATA part case
% 

%% S2P - Serial to Parallel conversion

S2P = zeros(N_SYM,N_SD); % DATA subcarriers initialized by zeros

k=0;

for i=1:N_SYM
    for j=1:N_SD
        S2P(i,j) = ofdm_input(k*N_SD+j);
    end
    k = k + 1;
end

%% Subcarriers initialization 

if ( strcmp(FORMAT,'non-HT') )
    
    % copy subcarriers for non-HT duplicate by N_SYM times 
    subcarriers = repmat(subcarriers_non_HT_dup,[N_SYM 1]);

    for n = 1:N_SYM

        subcarriers(n,pilot1_indexes) = P(1) * pn(n);
        subcarriers(n,pilot2_indexes) = P(2) * pn(n);
        subcarriers(n,pilot3_indexes) = P(3) * pn(n);
        subcarriers(n,pilot4_indexes) = P(4) * pn(n);

    end
    
elseif( strcmp(FORMAT,'VHT') )
    
    % copy subcarriers for VHT by N_SYM times 
    subcarriers = repmat(subcarriers_VHT,[N_SYM 1]);
    
    for n = 1:N_SYM

        subcarriers(n,pilot1_indexes) =  pilots(n,1);
        subcarriers(n,pilot2_indexes) = pilots(n,2);
        subcarriers(n,pilot3_indexes) = pilots(n,3);
        subcarriers(n,pilot4_indexes) = pilots(n,4);
        subcarriers(n,pilot5_indexes) = pilots(n,5);
        subcarriers(n,pilot6_indexes) = pilots(n,6);
        
    end
    
end

% data mapping

for n = 1:N_SYM
    j=1; % index of S2P
    for i = 1:IFFT_size
        if ( sum(i ~= pilot_indexes) == N_SP ) 
            % if subcarrier index is not pilot
            if (subcarriers(n,i) ~= 0)
                % if subcarrier is not NULL or DC
                subcarriers(n,i) = S2P(n,j);
                j = j + 1; % increment index of S2P
            end
        end
    end
end

subcarriers = ifftshift(subcarriers,2); % shifted IFFT block mapping

%% CP insertion

CP_length = T_GI;
CP_start = BCU_W - CP_length + 1; % start index
CP_stop = BCU_W; % stop index

% apply IFFT operation on subcarriers
for n = 1:N_SYM
    subcarriers(n,:) = ifft(subcarriers(n,:),IFFT_size);
end

% fit subcarriers into one BCU (adding effective guard band)
subcarriers_per_BCU = zeros(N_SYM,BCU_W);
guard_band_nulls = BCU_W - IFFT_size;

for n = 1:N_SYM
    subcarriers_per_BCU(n,:) = [zeros(1,guard_band_nulls/2) subcarriers(n,:) zeros(1,guard_band_nulls/2)];
end

% cycle prefix insertion
subcarriers_per_BCU_with_prefix = zeros(N_SYM,( BCU_W + CP_length ));

for n=1:N_SYM
    subcarriers_per_BCU_with_prefix(n,:) = [subcarriers_per_BCU(n,CP_start:CP_stop) subcarriers_per_BCU(n,:)];
end

col = BCU_W + T_GI; % number of columns
row = N_SYM; % number of rows

if ( strcmp(FORMAT,'VHT') )
    extra_zeros_flag = 0; % flag is set 0 because of not extra nulls necessary in VHT
end

%% P2S - converts full OFDM symbols to serial OFDM sybols stream 
ofdm_output=[];

for i=1:row
    ofdm_output = [ofdm_output subcarriers_per_BCU_with_prefix(i,:)];
end
   
end