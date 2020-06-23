function [ descrambler_output ] = descrambler( descrambler_input, s, half_byte)
%Function descrambles data input and returns descrambled data output
%descrambler_input - [bits] data input array
%s - descrambler initialization, value could be max. 127 "1111111" 
%descrambler_output - [bits] data output array
%s=127; % initialization sequence of descrambler (all 7 bits "1111111")

if ( size( half_byte , 2 ) == 2 ) % case if half byte size is equal to 2
    if ( (half_byte(end) > 1 ) || ( half_byte(1) > 1 ) ) % case if half byte is zeros
        extra_zeros = zeros(1, ( half_byte(end) - half_byte(1) ) + 1);
        descrambler_input = [descrambler_input( 1 : ( half_byte(1) - 1 ) ) extra_zeros descrambler_input( half_byte(1) : end ) ];
    else % case if half byte is not zeros
        descrambler_input = [descrambler_input half_byte];
    end
elseif ( size( half_byte , 2 ) == 4 ) % case if half byte is not zeros
    descrambler_input = [descrambler_input half_byte];
elseif ( size( half_byte , 2 ) == 6 ) % case if half byte is not zeros
    descrambler_input = [descrambler_input half_byte];
end

%% Data bits conversion to data bytes
j=1;
bytes=zeros(length(descrambler_input)/8,8);
for j=1:(length(descrambler_input)/8)
    for i=1:8
        bytes(j,i) = descrambler_input(i+(j-1)*8);
    end
    j=j+1;
end
bytes = bi2de(bytes,'left-msb');
descrambler_input = bytes';

%% Descramble data bytes

idata=zeros(size(descrambler_input));      % data field [B]

for k=1:size(descrambler_input,2);
    for i=1:8
        msb=bitxor(bitget(s,4),bitget(s,7)); %MSB bit calculation
        b4=bitget(s,4); %get 4th bit
        s=bitshift(s,-1); % shift scrambler register
        s=bitset(s,5,b4); % set 5th bit
        s=bitset(s,1,msb);% set MSB bit
        r=bitxor(bitget(descrambler_input(k),9-i),msb); % scrambled i-th data bit of data byte
        idata(k)=bitset(idata(k),9-i,r); % fill data field by scrambled bits
    end
end

%% Data bytes conversion to data bits

descrambler_output = [];

    bdata = de2bi(idata,8,'right-msb'); %  length(idata) x 8
    for i=length(idata):-1:1
        descrambler_output = [bdata(i,:), descrambler_output];
    end
    
    
    %if ( half_byte > -1 )
    %    descrambler_output = [descrambler_output(1:(end-8)) descrambler_output((end-3):end)];  
    %end
    
    if ( size( half_byte , 2 ) == 2 ) % case if half byte size is equal to 2
        if ( (half_byte(end) > 1 ) || ( half_byte(1) > 1 ) ) % case if half byte is zeros
            extra_zeros_length =  ( half_byte(end) - half_byte(1) ) + 1;
            if ( extra_zeros_length == 2 )
                descrambler_output = [descrambler_output(1:(end-8)) descrambler_output((end-5):end)]; 
            elseif ( extra_zeros_length == 4 )
                descrambler_output = [descrambler_output(1:(end-8)) descrambler_output((end-3):end)]; 
            elseif ( extra_zeros_length == 6 )
                descrambler_output = [descrambler_output(1:(end-8)) descrambler_output((end-1):end)]; 
            end         
        else % case if half byte is not zeros
            descrambler_output = [descrambler_output(1:(end-8)) descrambler_output((end-5):end)];
        end
    elseif ( size( half_byte , 2 ) == 4 ) % case if half byte is not zeros
        descrambler_output = [descrambler_output(1:(end-8)) descrambler_output((end-3):end)]; 
    elseif ( size( half_byte , 2 ) == 6 ) % case if half byte is not zeros
        descrambler_output = [descrambler_output(1:(end-8)) descrambler_output((end-1):end)]; 
    end

    
end
