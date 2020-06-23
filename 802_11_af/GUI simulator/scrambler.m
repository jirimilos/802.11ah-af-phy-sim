function [ scrambler_output, half_byte ] = scrambler( scrambler_input, s)
%Function scrambles data input and returns scrambled data output
%scrambler_input - [B|bits] data input array
%s - scrambler initialization, value could be max. 127 "1111111" 
%scrambler_output - [B] data output array
%s=127; % initialization sequence of scrambler (all 7 bits "1111111")

%% scrambler input detection 
flag = 0; % flag 0 -> scrambler_input are bytes
counter = 0; % counter of zeros and ones
for i = 1:length(scrambler_input)
    if ( scrambler_input(i) == 0 || scrambler_input(i) == 1)
        counter = counter + 1;
    end
end

if ( counter == length(scrambler_input) )
    flag = 1; % flag 1 -> scrambler_input are bits
end

if (flag == 1)
    
    % add half byte zeros without changing data for scrambler processing
    if ( mod(length(scrambler_input),8) == 8*0.5)
        extra_zeros = zeros(1,8*0.5);
        scrambler_input = [scrambler_input(1:(end-4)) extra_zeros scrambler_input((end-3):end)];
        half_byte_flag = 0;
    elseif ( mod(length(scrambler_input),8) == 8*0.25)
        extra_zeros = zeros(1,8*0.75);
        scrambler_input = [scrambler_input(1:(end-6)) extra_zeros scrambler_input((end-5):end)];
        half_byte_flag = 1;
    elseif ( mod(length(scrambler_input),8) == 8*0.75)
        extra_zeros = zeros(1,8*0.25);
        scrambler_input = [scrambler_input(1:(end-2)) extra_zeros scrambler_input((end-1):end)];
        half_byte_flag = 2;
    else
        half_byte_flag = -1;
        half_byte = half_byte_flag;
    end
    
    j=1;
    bytes=zeros(length(scrambler_input)/8,8);
    for j=1:(length(scrambler_input)/8)
        for i=1:8
            bytes(j,i) = scrambler_input(i+(j-1)*8);
        end
        j=j+1;
    end
    bytes = bi2de(bytes);
    scrambler_input = bytes';
% else
%     scrambler_input = scrambler_input;
end

idata=zeros(size(scrambler_input));      % data field [B]
for k=1:size(scrambler_input,2);
    for i=1:8
        msb=bitxor(bitget(s,4),bitget(s,7)); %MSB bit calculation
        b4=bitget(s,4); %get 4th bit
        s=bitshift(s,-1); % shift scrambler register
        s=bitset(s,5,b4); % set 5th bit
        s=bitset(s,1,msb);% set MSB bit
        r=bitxor(bitget(scrambler_input(k),9-i),msb); % scrambled i-th data bit of data byte
        idata(k)=bitset(idata(k),9-i,r); % fill data field by scrambled bits
    end
end
%scrambler_output = idata; % scrambled data input [B]

%% Data bytes conversion to data bits

scrambler_output = [];

    bdata = de2bi(idata,8,'left-msb'); %  length(idata) x 8
    for i=length(idata):-1:1
        scrambler_output = [bdata(i,:), scrambler_output];
    end
    
    if ( half_byte_flag == 0) % case for removing extra zeros and recording position
        counter = 0;
        start_index = 0; end_index = 0;
        half_byte = [start_index end_index];
        zeros_indexes = find( scrambler_output == 0 );
        for i = 2:length(zeros_indexes)
            if ( ( zeros_indexes(i) - zeros_indexes(i-1) ) == 1 )
                counter = counter + 1;
            else
                counter = 0;
            end
            if ( counter == 3)
                half_byte(1) = zeros_indexes(i-3);
                half_byte(2) = zeros_indexes(i);
                break;
            end
        end
    elseif ( half_byte_flag == 2)
        counter = 0;
        start_index = 0; end_index = 0;
        half_byte = [start_index end_index];
        zeros_indexes = find( scrambler_output == 0 );
        for i = 2:length(zeros_indexes)
            if ( ( zeros_indexes(i) - zeros_indexes(i-1) ) == 1 )
                counter = counter + 1;
            else
                counter = 0;
            end
            if ( counter == 1)
                half_byte(1) = zeros_indexes(i-1);
                half_byte(2) = zeros_indexes(i);
                break;
            end
        end
    elseif ( half_byte_flag == 1)
        counter = 0;
        start_index = 0; end_index = 0;
        half_byte = [start_index end_index];
        zeros_indexes = find( scrambler_output == 0 );
        for i = 2:length(zeros_indexes)
            if ( ( zeros_indexes(i) - zeros_indexes(i-1) ) == 1 )
                counter = counter + 1;
            else
                counter = 0;
            end
            if ( counter == 5)
                half_byte(1) = zeros_indexes(i-5);
                half_byte(2) = zeros_indexes(i);
                break;
            end
        end
        
    end
    
    if (half_byte_flag > -1 )
        if ( half_byte == [0 0] )
            warning('Half byte zeros not found!');
            if ( half_byte_flag == 0)
                half_byte = scrambler_output( ( end - 3 ) : end ); % gets the last part of scrambler output
                scrambler_output = scrambler_output( 1 : ( end - 4 ) ); % remove the last part of scrambler output
            elseif ( half_byte_flag == 2)
                half_byte = scrambler_output( ( end - 1 ) : end ); % gets the last part of scrambler output
                scrambler_output = scrambler_output( 1 : ( end - 2 ) ); % remove the last part of scrambler output
            elseif ( half_byte_flag == 1)    
                half_byte = scrambler_output( ( end - 5 ) : end ); % gets the last part of scrambler output
                scrambler_output = scrambler_output( 1 : ( end - 6 ) ); % remove the last part of scrambler output
            end
        else
            scrambler_output = [scrambler_output( 1 : ( half_byte(1) - 1 ) ) scrambler_output( ( half_byte(2) + 1 ) : end ) ];
        end
    end
    
end

