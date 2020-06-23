function [ bytes_output ] = bits2bytes( bits_input )
% Function converts binary input to bytes
% bits_input [bits] - binary input
% bytes_output [B] - bytes output

j=1;
bytes=zeros(length(bits_input)/8,8);
for j=1:(length(bits_input)/8)
    for i=1:8
        bytes(j,i) = bits_input(i+(j-1)*8);
    end
    j=j+1;
end
bytes = bi2de(bytes);
bytes_output = bytes';

end

