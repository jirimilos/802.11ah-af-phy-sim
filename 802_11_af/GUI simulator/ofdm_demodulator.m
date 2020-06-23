function [ ofdm_output ] = ofdm_demodulator( ofdm_input, guard_interval)
%Function process OFDM demodulation (removes GI interval) for OFDM symbols
%and returns complex subcarriers output
%ofdm_input - [complex numbers] data input to OFDM demodulator
%ofdm_output - [complex numbers] data output from OFDM demodulator

N_ST = 52; % subcarriers in total
FFT_size = 64; % size of IFFT window



switch(guard_interval)
    case 'GI'
        T_GI = floor(FFT_size/4); % GI duration
    case 'GIS' 
        T_GI = floor(FFT_size/8); % GI short duration
    otherwise error('Incorrect guard interval! Could be GI/GIS')
end

%% S2P - converts serial OFDM symbols stream to parallel complete OFDM symbols  
ofdm_symbol_length=FFT_size + T_GI; % length of 1 complete OFDM symbol
ofdm_symbols_count = length(ofdm_input) / ofdm_symbol_length; % count of OFDM symbols in OFDM signal
ofdm_symbols_full = zeros(ofdm_symbols_count,ofdm_symbol_length); % matrix of complete OFDM symbols
k=0;
for j=1:ofdm_symbols_count
    for i=1:ofdm_symbol_length
        ofdm_symbols_full(j,i) = ofdm_input(ofdm_symbol_length*k+i);
    end
    k = k + 1;
end

%% Remove GI

start_ofdm_symbol = T_GI + 1; % start position for OFDM symbol
ofdm_symbols = ofdm_symbols_full(:,start_ofdm_symbol:end); % OFDM symbols without GI


%% FFT 
fft_output = zeros(ofdm_symbols_count,FFT_size); % matrix of complete OFDM symbols
for i=1:ofdm_symbols_count
    fft_output(i,:)=fft(ofdm_symbols(i,:)); % OFDM symbol FFT transformation
end

%% Remove NULL and DC subcarriers
fft_subcarriers = fft_output';
%fft_subcarriers = round(real(fft_output))'; % define all subcarriers real part + round numbers
DPS = ofdm_symbols_count; % number of data per symbol
subcarriers = zeros(N_ST,DPS);

    for i=1:FFT_size
        if (i > 1) && (i < 28) % subcarrier indexes from 1 to 26
                                   % get from FFT block inputs from 1 to 26
            subcarriers(i+25,:) = fft_subcarriers(i,:);                     
        elseif (i > 38) % subcarrier indexes from -26 to -1 
                                % get from FFT block inputs from 38 to 63 
            subcarriers(i-38,:) = fft_subcarriers(i,:); 
        end
    end
    
%% Remove PILOTS subcarriers
modulated_data = subcarriers([1:5,7:19,21:32,34:46,48:end],:); % subcarriers 48*DATA modulated

%% P2S conversion
data_length = size(modulated_data,1) * size(modulated_data,2); % length of modulated DATA

ofdm_output=[];
N_ST = size(modulated_data,1);

k=0;

for i=1:DPS % convert parallel stream to serial
    for j=1:N_ST
        ofdm_output = [ofdm_output modulated_data(j,i)];
    end
    k = k + 1;
end

