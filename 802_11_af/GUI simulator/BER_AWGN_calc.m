function [ BER_output ] = BER_AWGN_calc( SNR_vector, data_rate )
%Functions calculates theoretical BER vector using SNR and data rate
% SNR_vector - Signal to Noise ratio vector
% data_rate - Data rate used in IEEE 802.11

    switch(data_rate)
        case {'1.1','1.6'} % BPSK
            mod_type = 'psk'; % type of modulation
            M = 2;  % number of symbols
            N_CBPS = log2(M); % number of bits per symbol
            % calulate EbNo vector
            EbNo_vector = SNR_vector - 10 * log10(N_CBPS); 
            BER_output = berawgn(EbNo_vector,mod_type,M,'nondiff');
        case {'2.1','3.2'} % QPSK
            mod_type = 'psk'; % type of modulation
            M = 4; % number of symbols
            N_CBPS = log2(M); % number of bits per symbol
            % calulate EbNo vector
            EbNo_vector = SNR_vector - 10 * log10(N_CBPS); 
            BER_output = berawgn(EbNo_vector,mod_type,M,'nondiff');
         case {'4.3','6.4'} % 16-QAM
            mod_type = 'qam'; % type of modulation
            M = 16; % number of symbols
            N_CBPS = log2(M); % number of bits per symbol
            % calulate EbNo vector
            EbNo_vector = SNR_vector - 10 * log10(N_CBPS);
            BER_output = berawgn(EbNo_vector,mod_type,M);
        case {'8.5','9.6'} % 64-QAM
            mod_type = 'qam'; % type of modulation
            M = 64; % number of symbols
            N_CBPS = log2(M); % number of bits per symbol
            % calulate EbNo vector
            EbNo_vector = SNR_vector - 10 * log10(N_CBPS);
            BER_output = berawgn(EbNo_vector,mod_type,M);
        otherwise
            error('Data rate is incorrect!')
    end
    
end

