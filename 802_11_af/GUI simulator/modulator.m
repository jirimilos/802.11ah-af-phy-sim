function [ modulator_output ] = modulator( modulator_input, modulation)
%Function modulates data interleaved data input bits to
%BPSK/QPSK/16-QAM/64-QAM/256-QAM
%modulator_input - [bits] data input bits from interleaver
%modulation - [string] 'BPSK'/'QPSK'/'16-QAM'/'64-QAM'/'256-QAM'
%modulator_output - modulated data output symbols

% Subcarrier modulation mapping
bits_input = reshape(modulator_input',1,size(modulator_input,1)*size(modulator_input,2)); % FEC coded and interleaved bit serial input before modulation

% BPSK, QPSK, 16-QAM, 64-QAM, 256-QAM modulations
switch(modulation)
    case 'BPSK' 
        K_MOD = 1; % normalization factor 1
        N_BPSC = 1; % 1 coded bits per subcarrier
%         BPSK=[];
%         
%         for i=1:length(bits_input)        
%             if ( bits_input(i) == 0)
%                 I=-1; Q=0;
%             else
%                 I=1; Q=0;
%             end
%             d=complex(I,Q)*K_MOD; % output symbol
%             BPSK=[BPSK,d]; % BPSK modulated bits input
%         end
        BPSKModulator = comm.BPSKModulator;
        BPSK_modulated = step(BPSKModulator,bits_input'); % BPSK modulation
        BPSK_real = (-1)*real(BPSK_modulated); % BPSK corection according to the standard IEEE 802.11
        BPSK = complex(BPSK_real,imag(BPSK_modulated)); % BPSK output
        modulator_output = zeros(1,length(BPSK));
        
        for i=1:length(BPSK)
            modulator_output(i) = BPSK(i,1); 
        end

    case 'QPSK' 
%         K_MOD = 1/sqrt(2); % normalization factor 1/sqrt(2)
%         N_BPSC = 2; % 2 coded bits per subcarrier
%         QPSK=[];
%         for i=1:2:length(bits_input)
%             switch(bits_input(i))
%                 case 0
%                     if ( bits_input(i+1) == 0)
%                         I=-1; Q=-1; % B0B1 -> 00
%                     else
%                         I=-1; Q=1; %B0B1 -> 01
%                     end
%                 case 1
%                     if ( bits_input(i+1) == 0)
%                         I=1; Q=-1; %B0B1 -> 10
%                     else
%                         I=1; Q=1; %B0B1 -> 11
%                     end
%                 otherwise error ('Incorrect bit value! Bit could be 0/1')
%             end
%     
%             d=complex(I,Q)*K_MOD; % output symbol
%             QPSK=[QPSK,d]; % QPSK modulated bits input
% 
%         end
%         modulator_output = QPSK;

        K_MOD = 1/sqrt(2); % normalization factor 1/sqrt(2)
        N_BPSC = 2; % 2 coded bits per subcarrier
        
        h = modem.qammod('M',4);
        h.SymbolOrder = 'user-defined';
        h.SymbolMapping = [1 3 0 2];
        h.InputType = 'bit';
        
        data_symbols = modulate(h, bits_input');
        data_symbols = reshape(data_symbols,1,size(data_symbols,1)*size(data_symbols,2));
        modulator_output = data_symbols.*K_MOD;
        
    case '16-QAM' 
        K_MOD = 1/sqrt(10); % normalization factor 1/sqrt(10)
        N_BPSC = 4; % 4 coded bits per subcarrier
        
        h = modem.qammod('M',16);
        h.SymbolOrder = 'user-defined';
        h.SymbolMapping = [2 3 1 0 6 7 5 4 14 15 13 12 10 11 9 8];
        h.InputType = 'bit';
        
        data_symbols = modulate(h, bits_input');
        data_symbols = reshape(data_symbols,1,size(data_symbols,1)*size(data_symbols,2));
        modulator_output = data_symbols.*K_MOD;
        
    case '64-QAM' 
        K_MOD = 1/sqrt(42); % normalization factor 1/sqrt(42)
        N_BPSC = 6; % 6 coded bits per subcarrier
        
        h = modem.qammod('M',64);
        h.SymbolOrder = 'user-defined';
        h.SymbolMapping = [4 5 7 6 2 3 1 0 12 13 15 14 10 11 9 8 28 29 31 30 26 27 25 24 20 21 23 22 18 19 17 16 52 53 55 54 50 51 49 48 60 61 63 62 58 59 57 56 44 45 47 46 42 43 41 40 36 37 39 38 34 35 33 32];
        h.InputType = 'bit';
        
        data_symbols = modulate(h, bits_input');
        data_symbols = reshape(data_symbols,1,size(data_symbols,1)*size(data_symbols,2));
        modulator_output = data_symbols.*K_MOD;
        %modulator_output = reshape(modulator_output,1,size(modulator_output,1)*size(modulator_output,2));
    
    case '256-QAM'
        K_MOD = 1/sqrt(170)
        N_BPSC = 8; % 8 coded bits per subcarrier
        
        h = modem.qammod('M',256);
        h.SymbolOrder = 'gray';
        h.InputType = 'bit';
        
        data_symbols = modulate(h, bits_input');
        data_symbols = reshape(data_symbols,1,size(data_symbols,1)*size(data_symbols,2));
        modulator_output = data_symbols.*K_MOD;
        
    otherwise error('Modulation incorrect! Could be BPSK/QPSK/16-QAM/64-QAM/256-QAM')
end

%figure
%h=scatterplot(BPSK)
%plot(real(BPSK),imag(BPSK),'o')
%grid on

end
