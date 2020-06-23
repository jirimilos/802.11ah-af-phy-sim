function [ demodulator_output ] = demodulator( demodulator_input, modulation)
%Function demodulates data with BPSK/QPSK/16-QAM/64-QAM/256-QAM modulation
%demodulator_input - [symbols] modulated data symbols
%modulation - [string] 'BPSK'/'QPSK'/'16-QAM'/'64-QAM'/'256-QAM'
%demodulator_output - [bits] data bits to deinterleaver

% Subcarrier modulation mapping

% BPSK, QPSK, 16-QAM, 64-QAM, 256-QAM modulations
switch(modulation)
    case 'BPSK' 
        K_MOD = 1; % normalization factor 1
        N_BPSC = 1; % 1 coded bits per subcarrier
%         data = []; % demodulator output
%         
%         for i=1:length(demodulator_input)
%             symbol = round(real(demodulator_input(i)))/K_MOD; % symbol without normalization factor
%             if (symbol == 1) % if symbol 1 -> bit = 1
%                 bit = 1;
%             else % if symbol -1 -> bit = 0
%                 bit = 0; 
%             end
%             data=[data,bit]; % demodulated output bits
%         end
%         demodulator_output = data;

        demodulator_input = demodulator_input./K_MOD; % symbols without normalization factor    

        BPSK_modulated = zeros(length(demodulator_input),1);
        
        for i=1:length(demodulator_input)
            BPSK_modulated(i,1) = demodulator_input(i); 
        end
        
        BPSK_real = (-1)*real(BPSK_modulated); % BPSK corection according to the standard IEEE 802.11
        BPSK = complex(BPSK_real,imag(BPSK_modulated)); % BPSK output
        
        
        BPSKDemodulator = comm.BPSKDemodulator;
        BPSK_demodulated = step(BPSKDemodulator,BPSK); % BPSK demodulation
        
        demodulator_output = [];
        for i=1:length(demodulator_input)
            demodulator_output(i) = BPSK_demodulated(i,1);
        end
        %demodulator_output = BPSK_demodulated;
        
    case 'QPSK' 
%         K_MOD = 1/sqrt(2); % normalization factor 1/sqrt(2)
%         N_BPSC = 2; % 2 coded bits per subcarrier
%         data=[]; % demodulator output
%         for k=1:length(demodulator_input)
%             symbol = round(demodulator_input(k)/K_MOD); %  symbol without normalization factor
%             
%             if symbol == complex(-1,-1)
%                 data = [data, [0 0]]; % 00 bits
%             elseif symbol == complex(-1,1)
%                 data = [data, [0 1]]; % 01 bits
%             elseif symbol == complex(1,-1)
%                 data = [data, [1 0]]; % 10 bits
%             elseif symbol == complex(1,1)
%                 data = [data, [1 1]]; % 11 bits
%             else
%                 error ('Incorrect symbol value!')
%             end
% 
%         end
%         demodulator_output = data;

        K_MOD = 1/sqrt(2); % normalization factor 1/sqrt(2)
        N_BPSC = 2; % 2 coded bits per subcarrier
        
        demodulator_input = demodulator_input./K_MOD; % symbols without normalization factor
        
        h = modem.qamdemod('M',4);
        h.SymbolOrder = 'user-defined';
        h.SymbolMapping = [1 3 0 2];
        h.OutputType = 'bit';
        
        data_sybols = [];
        for i=1:length(demodulator_input)
            data_symbols(i,1) = demodulator_input(i);
        end
        
        data_bits = demodulate(h, data_symbols);
        demodulator_output = data_bits';

    case '16-QAM' 
        K_MOD = 1/sqrt(10); % normalization factor 1/sqrt(10)
        N_BPSC = 4; % 4 coded bits per subcarrier
        
        demodulator_input = demodulator_input./K_MOD; % symbols without normalization factor
        
        h = modem.qamdemod('M',16);
        h.SymbolOrder = 'user-defined';
        h.SymbolMapping = [2 3 1 0 6 7 5 4 14 15 13 12 10 11 9 8];
        h.OutputType = 'bit';
        
        data_sybols = [];
        for i=1:length(demodulator_input)
            data_symbols(i,1) = demodulator_input(i);
        end
        
        data_bits = demodulate(h, data_symbols);
        demodulator_output = data_bits';
        
    case '64-QAM' 
        K_MOD = 1/sqrt(42); % normalization factor 1/sqrt(42)
        N_BPSC = 6; % 6 coded bits per subcarrier
       
        demodulator_input = demodulator_input./K_MOD; % symbols without normalization factor
        
        h = modem.qamdemod('M',64);
        h.SymbolOrder = 'user-defined';
        h.SymbolMapping = [4 5 7 6 2 3 1 0 12 13 15 14 10 11 9 8 28 29 31 30 26 27 25 24 20 21 23 22 18 19 17 16 52 53 55 54 50 51 49 48 60 61 63 62 58 59 57 56 44 45 47 46 42 43 41 40 36 37 39 38 34 35 33 32];
        h.OutputType = 'bit';
        
        data_sybols = [];
        for i=1:length(demodulator_input)
            data_symbols(i,1) = demodulator_input(i);
        end
        
        data_bits = demodulate(h, data_symbols);
        demodulator_output = data_bits';
        
    case '256-QAM'
        K_MOD = 1/sqrt(170)
        N_BPSC = 8; % 8 coded bits per subcarrier
        
        demodulator_input = demodulator_input./K_MOD; % symbols without normalization factor
        
        h = modem.qamdemod('M',256);
        h.SymbolOrder = 'gray';
        h.OutputType = 'bit';
        
        data_sybols = [];
        for i=1:length(demodulator_input)
            data_symbols(i,1) = demodulator_input(i);
        end
        
        data_bits = demodulate(h, data_symbols);
        demodulator_output = data_bits';
        
    otherwise error('Modulation incorrect! Could be BPSK/QPSK/16-QAM/64-QAM/256-QAM')
end

%figure
%h=scatterplot(BPSK)
%plot(real(BPSK),imag(BPSK),'o')
%grid on

end
