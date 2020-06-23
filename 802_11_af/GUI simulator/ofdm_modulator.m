function [ ofdm_output, ifft_output ] = ofdm_modulator( ofdm_input, guard_interval, cyclic_prefix)
%Function process OFDM modulation of complex modulated subcarriers
%and returns OFDM complex modulated output
%ofdm_input - [complex numbers] data input to OFDM modulator
%ofdm_output - [complex numbers] data output from OFDM modulator

N_ST = 52; % subcarriers in total
IFFT_size = 64; % size of IFFT window
N_SYM = length(ofdm_input)/48; % number of OFDM symbols
% pn sequence of pilots polarity for 127 OFDM symbols
p127 = [1,1,1,1, -1,-1,-1,1, -1,-1,-1,-1, 1,1,-1,1, -1,-1,1,1, -1,1,1,-1, 1,1,1,1, 1,1,-1,1, 1,1,-1,1, 1,-1,-1,1, 1,1,-1,1, -1,-1,-1,1, -1,1,-1,-1, 1,-1,-1,1, 1,1,1,1, -1,-1,1,1, -1,-1,1,-1, 1,-1,1,1, -1,-1,-1,1, 1,-1,-1,-1, -1,1,-1,-1, 1,-1,1,1, 1,1,-1,1, -1,1,-1,1, -1,-1,-1,-1, -1,1,-1,1, 1,-1,1,-1, 1,1,1,-1, -1,1,-1,-1, -1,1,1,1, -1,-1,-1,-1, -1,-1,-1];
% pilots sequence of pilots -7, -21, +7, +21
P = [1 1 1 -1]; 
p127_coeff = ceil(N_SYM/127);
if ( p127_coeff == 1)
    pn = p127;
else % extend pn sequency according to number of OFDM symbols
    pn = repmat(p127,1,p127_coeff);
end

switch(cyclic_prefix)
    case 0
        CP = 0; % cyclic prefix disable
    case 1 
        CP = 1; % cyclic prefix enable
    otherwise error('Incorrect cyclic prefix! Could be 0/1(disable/enable)')
end

switch(guard_interval)
    case 'GI'
        T_GI = floor(IFFT_size/4); % GI duration
    case 'GIS' 
        T_GI = floor(IFFT_size/8); % GI short duration
    otherwise error('Incorrect guard interval! Could be GI/GIS')
end

%% S2P - Serial to Parallel conversion

data_length = length(ofdm_input); % length of modulated DATA symbols (complex numbers)
N_SD = 48; % number of streams represented by DATA symbols (complex numbers)
DPS= data_length/N_SD; % length of DATA symbols per 1 subcarrier
S2P = zeros(N_SD,DPS); % 48 DATA subcarriers initialized by zeros 

k=0;
for i=1:DPS % convert serial stream to parallel
    for j=1:N_SD
        S2P(j,i) = ofdm_input(k*N_SD+j);
    end
    k = k + 1;
end

%subcarriers = S2P'; % 48 DATA subcarriers

%% Subcarriers initialization 
% Subcarriers (48 x DATA + 4 x PILOT) 
% size allocation
subcarriers = zeros(N_ST,DPS); 

% Subcarriers (48 x DATA + 4 x PILOT) 
% initialization 

% pilot initialization
for n=1:DPS
    
    subcarriers(6,n) = P(1) * pn(n);
    subcarriers(20,n) = P(2) * pn(n);
    subcarriers(33,n) = P(3) * pn(n);
    subcarriers(47,n) = P(4) * pn(n);
    
end

j=1; % index of DATA subcarrier

% data initialization
for i=1:N_ST
    if ( ( i ~= 6 ) && ( i ~= 20 ) && ( i ~= 33 ) && ( i ~=47 ) )
        % Data with indexes from -26 to +26 except 4 x pilots
        subcarriers(i,:) = S2P(j,:);
        j = j + 1;
    end
end

% % data initialization
% for i=1:N_ST
%     if i==6
%         % Pilot with index -21
%         subcarriers(i,:) = ones(1,DPS); % 1
%     elseif i == 20
%         % Pilot with index -7
%         subcarriers(i,:) = ones(1,DPS); % 1
%     elseif i == 33
%         % Pilot with index +7
%         subcarriers(i,:) = ones(1,DPS); % 1
%     elseif i == 47
%         % Pilot with index +21
%         subcarriers(i,:) = (-1)*ones(1,DPS); % -1
%     else
%         % Data with indexes from -26 to +26 except 4 x pilots
%         subcarriers(i,:) = S2P(j,:);
%         j = j + 1;
%     end
% end

%% IFFT block mapping
% for 48 x data + 4 x pilot + 1 x DC + 11 x null subcarriers
    ifft_subcarriers = zeros(IFFT_size,DPS); % 64 * DPS
    
    for i=1:N_ST
        if (i > 0) && (i < 27) % subcarrier indexes from -26 to -1
                                   % go to IFFT block inputs from 38 to 63
            ifft_subcarriers(i+38,:) = subcarriers(i,:);                     
        elseif (i > 26) && (i < 53) % subcarrier indexes from 1 to 26 
                                % go to IFFT block inputs from 1 to 26 
            ifft_subcarriers(i-25,:) = subcarriers(i,:);
        end
    end
    

%% GI/CP insertion

% (1) situation when cyclic prefix enabled
if CP == 1
    
    CP_len = T_GI;
    CP_start = IFFT_size - CP_len + 1; % start index
    CP_stop = IFFT_size; % stop index
    
    % allocated IFFT block Inverse Fourier Transformation 
    ifft_data = zeros(DPS, IFFT_size); % empty matrix
    for i=1:DPS
        ifft_data(i,:) = [ifft(ifft_subcarriers(:,i)')];
    end   
    
    % CP insertion
    for i=1:DPS
        ifft_output(i,:) = [ifft_data(i,CP_start:CP_stop) ifft_data(i,:)];
    end

% (2) situation when cyclic prefix disabled
else 
    GI = zeros(1,T_GI); % guard interval block
    
    for i=1:DPS
        ifft_output(i,:) = [GI ifft(ifft_subcarriers(:,i)')];
    end

end

col = IFFT_size + T_GI; % number of columns
row = DPS; % number of rows


%% P2S - converts full OFDM symbols to serial OFDM sybols stream 
ofdm_output=[];

for i=1:row
    ofdm_output = [ofdm_output ifft_output(i,:)];
end
    ofdm_signal_length = col * row; % full length of OFDM signal
%    ofdm_output = reshape(ifft_output, 1, ofdm_signal_length); % OFDM signal to transmit
