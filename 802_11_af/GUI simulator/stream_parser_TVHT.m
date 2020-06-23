function [ output_streams ] = stream_parser_TVHT(coded_data)

N_BPSCS = 1; % number of coded bits per symbol per spatial stream
N_SS=1; % number of spatial streams
N_ES=1; % number of encoders
N_CBPS = 216; % number of coded bits per symbol
N_SD = 108;
N_CBPSS = N_SD; % number of coded bit per symbol per spatial stream

s = max(1,(N_BPSCS/2)); % number of bits assigned to a single axis of a constellation point in spatial stream
S = N_SS * s; % sum of theese over all streams

N_Block = floor(N_CBPS/(N_ES*S));
M = (N_CBPS - ( N_Block * N_ES * S) ) / ( s * N_ES); 


%k_apost = k - N_Block * N_ES * s;
%L = floor(k_apost/s) * N_SS + (i_ss - 1);

k_max1 = N_Block * N_ES * s - 1;
k2 = N_Block * N_ES * s;
k_max2 = N_CBPSS-1; 


j = [];
i = [];

% calculate j
for k = 0:k_max1
    j_act = mod(floor(k/s),N_ES);
    j = [j j_act];
end
for k = k2:k_max2
    j_act = floor(L/M);
    j = [j j_act];
end


% calculate i
for i_ss=1:N_SS
    for k=0:k_max1
        i_act = (i_ss - 1) * s + S * floor(k/(N_ES*s)) + mod(k,s);
        i = [i i_act];
    end
end

for i_ss=1:N_SS
    for k=k2:k_max2
        k_apost = k - N_Block * N_ES * s;
        L = floor(k_apost/s) * N_SS + (i_ss - 1);
        i_act = mod(L,M) * s + N_Block * S + mod(k,s);
        i = [i i_act];
    end
end

end

