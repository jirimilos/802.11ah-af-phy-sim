function [ interleaver_output1, interleaver_output2 ] = interleaver( interleaver_input, CH_BANDWIDTH, transmission_type, modulation)
%Function interleaves data coded input and returns interleaved data output
%interleaver_input - [bits] data FEC coded input
%modulation - [string] BPSK/QPSK/16-QAM/64-QAM modulation
%interleaver_output1 - [bits] data output after 1st permutation
%interleaver_output2 - [bits] data output after 2nd permutation

switch(transmission_type)
    case 'non-HT'
        switch(modulation)
            case 'BPSK' 
                N_CBPS = 48; % 48 coded bits per OFDM symbol
                N_BPSC = 1; % 1 coded bits per subcarrier
            case 'QPSK' 
                N_CBPS = 96; % 96 coded bits per OFDM symbol
                N_BPSC = 2; % 2 coded bits per subcarrier
            case '16-QAM' 
                N_CBPS = 192; % 192 coded bits per OFDM symbol
                N_BPSC = 4; % 4 coded bits per subcarrier
            case '64-QAM' 
                N_CBPS = 288; % 288 coded bits per OFDM symbol
                N_BPSC = 6; % 6 coded bits per subcarrier
            otherwise error('Modulation incorrect! Could be BPSK/QPSK/16-QAM/64-QAM')
        end
    case 'VHT'
        switch(modulation)
            case 'BPSK' 
                N_CBPS = 108; % 108 coded bits per OFDM symbol
                N_BPSC = 1; % 1 coded bits per subcarrier
            case 'QPSK' 
                N_CBPS = 216; % 216 coded bits per OFDM symbol
                N_BPSC = 2; % 2 coded bits per subcarrier
            case '16-QAM' 
                N_CBPS = 432; % 432 coded bits per OFDM symbol
                N_BPSC = 4; % 4 coded bits per subcarrier
            case '64-QAM' 
                N_CBPS = 648; % 648 coded bits per OFDM symbol
                N_BPSC = 6; % 6 coded bits per subcarrier
            case '256-QAM' 
                N_CBPS = 864; % 864 coded bits per OFDM symbol
                N_BPSC = 8; % 8 coded bits per subcarrier
            otherwise error('Modulation incorrect! Could be BPSK/QPSK/16-QAM/64-QAM/256-QAM')
        end
end
        
switch(CH_BANDWIDTH)
    case 'TVHT_W' 
        N_COL = 18; % number of bits in interleaver block columns
        N_ROW = 6 * N_BPSC; % number of bits in interleaver row
        N_SEG = 1; % number of channel segments 
    case {'TVHT_2W','TVHT_W+W'} 
        N_CBPS = 2 * N_CBPS;
        N_COL = 27; % number of bits in interleaver block columns
        N_ROW = 8 * N_BPSC; % number of bits in interleaver row
        N_SEG = 2; % number of channel segments 
    case {'TVHT_4W','TVHT_2W+2W'}
        N_CBPS = 4 * N_CBPS;
        N_COL = 48; % number of bits in interleaver block columns
        N_ROW = 9 * N_BPSC; % number of bits in interleaver row
        N_SEG = 4; % number of channel segments 
    otherwise error('Incorrect channel bandwidth used! Could be TVHT_(W/2W/4W/W+W/2W+2W)')
end

s=max(N_BPSC/2,1); % s value 

switch( transmission_type)
    case 'VHT'
        N_SD = 108; % only for TVHT MODE 1 (N_SS=1)
        N_SD = N_SD * N_SEG; % number of complex data per BCU
        N_CBPSSI = N_CBPS; % number of coded bits per symbol per spatial stream per BCC interleaver block
        %N_BPSCS = N_BPSC; % number of coded bits per subcarrier per spatial stream for user X = 0
        
        %N_COL = 18; % number of bits in interleaver block columns
        %N_ROW = 6 * N_BPSCS; % number of bits in interleaver row
        
        %blocks=reshape(interleaver_input,length(interleaver_input)/N_CBPS,N_CBPS);
        blocks=reshape(interleaver_input,N_CBPS,length(interleaver_input)/N_CBPS)';
        number_of_blocks = length(interleaver_input)/N_CBPS; % number of input blocks
        
        index_i = [];
        II_1=blocks;
        for n=1:number_of_blocks
            for k=0:(N_CBPSSI-1)
                i = N_ROW * mod(k,N_COL) + floor(k/N_COL);
                II_1(n,k+1) = blocks(n,i+1);
                index_i = [index_i,i];
            end
        end
        
        index_j = [];
        II_2=II_1;
        for n=1:number_of_blocks
            for i=0:(N_CBPSSI-1)
                a = s*floor(i/s);
                b=floor((N_COL*i)/N_CBPSSI);
                j=a+mod((i+N_CBPSSI-b),s); % index after the second permutation
                II_2(n,i+1)=II_1(n,j+1);
                index_j = [index_j,j];
            end
        end        
        
        %index_j = reshape(index_j,N_CBPS,length(interleaver_input)/N_CBPS)';

        interleaver_output1 = II_1;
        interleaver_output2 = II_2;
    
    case 'non-HT'     

        II=reshape(interleaver_input,N_CBPS,length(interleaver_input)/N_CBPS)';

        number_of_blocks = length(interleaver_input)/N_CBPS; % number of input blocks
        index_i = [];
        II_1=II;
        for n=1:number_of_blocks
            for k=0:(N_CBPS-1)
                i=(N_CBPS/16)*mod(k,16)+floor(k/16);
                II_1(n,k+1)=II(n,i+1);
                index_i = [index_i,i];
            end
        end

        %index_i = reshape(index_i,N_CBPS,length(interleaver_input)/N_CBPS)';

        index_j = [];
        II_2=II_1;
        for n=1:number_of_blocks
            for i=0:(N_CBPS-1)
                a = s*floor(i/s);
                b=floor((16*i)/N_CBPS);
                j=a+mod((i+N_CBPS-b),s); % index after the second permutation
                %j=s*floor(i/s) + mod(i + N_CBPS - floor((16*i)/N_CBPS),s); % index after the first permutation
                II_2(n,i+1)=II_1(n,j+1);
                index_j = [index_j,j];
            end
        end
        index_j = reshape(index_j,N_CBPS,length(interleaver_input)/N_CBPS)';

        interleaver_output1 = II_1;
        interleaver_output2 = II_2;
end

end

