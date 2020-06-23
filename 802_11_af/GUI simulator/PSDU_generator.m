function [ data_output ] = PSDU_generator( N )
%Function generates random data bytes
%N - length of random data bytes
%data_output - output randomized data bytes with length N
data_output=randi(255,1,N); % generates random data bytes of 1 x N size
data_output=de2bi(data_output,8); % converts output to binary output matrix
data_length=size(data_output,1)*size(data_output,2); % data length
data_output=reshape(data_output',1,data_length); % PSDU data output
end