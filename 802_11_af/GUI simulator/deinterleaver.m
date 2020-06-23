function [ deinterleaver_output ] = deinterleaver( deinterleaver_input, CH_BANDWIDTH, transmission_type, modulation)
%Function deinterleaves data input and returns FEC coded output
%deinterleaver_input - [bits] matrix data_bits x N_SD data interleaved input
%modulation - [string] BPSK/QPSK/16-QAM/64-QAM modulation
%deinterleaver_output - [bits] FEC coded data output

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
           
        % reshape deinterleaver input to matrix with blocks length N_CBPSSI
        row=length(deinterleaver_input)/N_CBPSSI;
        col=N_CBPSSI;
        idata = zeros(row,col);
        k=0;
        for j=1:row
            for i=1:col
                idata(j,i) = deinterleaver_input(k*N_CBPSSI+i);
            end
            k=k+1;
        end

        number_of_blocks = size(idata,1); % number of input blocks
        %number_of_blocks = length(deinterleaver_input)/N_CBPS; % number of input blocks
        index_i = [];
        II_2=idata; % data bits matrix after 2nd permutation of interleaver
        II_1=zeros(number_of_blocks,N_CBPSSI);
        j=0;

        for n=1:number_of_blocks
            for j=0:(N_CBPSSI-1)
                a=j + floor((N_COL*j)/N_CBPSSI);
                b=s*floor(j/s);
                i=b+mod(a,s); % index i after the first permutation
                II_1(n,j+1)=II_2(n,i+1);
                index_i = [index_i,i];
            end
        end

        index_k = [];
        II=zeros(number_of_blocks,N_CBPSSI);
        for n=1:number_of_blocks
            for i=0:(N_CBPSSI-1)
                a=floor(i/N_ROW);
                b=N_COL*i;
                k=b-(N_CBPSSI-1)*a; % index k after the second permutation
                II(n,i+1)=II_1(n,k+1);
                index_k = [index_k,k];
            end
        end

        deinterleaver_output = reshape(II',1,length(index_k)); % output FEC coded data        
    case 'non-HT'
        % reshape deinterleaver input to matrix with blocks length N_CBPS
        row=length(deinterleaver_input)/N_CBPS;
        col=N_CBPS;
        idata = zeros(row,col);
        k=0;
        for j=1:row
            for i=1:col
                idata(j,i) = deinterleaver_input(k*N_CBPS+i);
            end
            k=k+1;
        end

        number_of_blocks = size(idata,1); % number of input blocks
        %number_of_blocks = length(deinterleaver_input)/N_CBPS; % number of input blocks
        index_i = [];
        II_2=idata; % data bits matrix after 2nd permutation of interleaver
        II_1=zeros(number_of_blocks,N_CBPS);
        j=0;

        for n=1:number_of_blocks
            for j=0:(N_CBPS-1)
                a=j + floor((16*j)/N_CBPS);
                b=s*floor(j/s);
                i=b+mod(a,s); % index i after the first permutation
                II_1(n,j+1)=II_2(n,i+1);
                index_i = [index_i,i];
            end
        end

        index_k = [];
        II=zeros(number_of_blocks,N_CBPS);
        for n=1:number_of_blocks
            for i=0:(N_CBPS-1)
                a = floor((16*i)/N_CBPS);
                b=16*i;
                k=b-(N_CBPS-1)*a; % index k after the second permutation
                II(n,i+1)=II_1(n,k+1);
                index_k = [index_k,k];
            end
        end

        deinterleaver_output = reshape(II',1,length(index_k)); % output FEC coded data
end

end