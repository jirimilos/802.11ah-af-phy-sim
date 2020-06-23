function [ coder_output ] = fec_coder( coder_input , rate)
%Function codes scrambled data input using convolution coder K=7
%coder_input - [B|bits] scramled data input array
%rate - [string] coder rate 1/2, 2/3, 3/4 or 5/6 for IEEE 802.11af
%coder_output - [bits] convolutionally coded data output [array]

switch(rate)
    case '1/2' 
        rate_vector = [1;1]; %R=1/2
    case '2/3' 
        rate_vector = [1;1;1;0]; %R=2/3
    case '3/4' 
        rate_vector = [1;1;1;0;0;1]; %R=3/4
    case '5/6' 
        rate_vector = [1;1;1;0;0;1;1;0;0;1]; %R=5/6
    otherwise error('Coding rate incorrect!')
end

K=7; % encoder length
G0=133; % generator polynomial G0
G1=171; % generator polynomial G1
t=poly2trellis(K,[G1 G0]); % define trellis

%% coder input detection 
flag = 0; % flag 0 -> coder_input are bytes
counter = 0; % counter of zeros and ones
for i = 1:length(coder_input)
    if ( coder_input(i) == 0 || coder_input(i) == 1)
        counter = counter + 1;
    end
end

if ( counter == length(coder_input) )
    flag = 1; % flag 1 -> coder_input are bits
end

binary_input=[];
if (flag == 0)
    bdata = de2bi(coder_input,8,'left-msb'); %  length(data_input) x 8
    for i=length(coder_input):-1:1
        binary_input = [bdata(i,:), binary_input];
    end
else
    binary_input = coder_input;
end

% if gpuDeviceCount == 0 % gpu device not available
    H = comm.ConvolutionalEncoder('TrellisStructure',t,...
        'TerminationMethod','Truncated',...
        'PuncturePatternSource','Property',...
        'PuncturePattern',rate_vector);
% else    
%     H = comm.gpu.ConvolutionalEncoder('TrellisStructure',t,...
%         'TerminationMethod','Truncated',...
%         'PuncturePatternSource','Property',...
%         'PuncturePattern',rate_vector,...
%         'NumFrames',1);
% end                              
                              
convolutional_encoder_output = step(H, binary_input');

coder_output = convolutional_encoder_output';


%coder_output = convenc(binary_input,t,rate_vector); % coded data input

end