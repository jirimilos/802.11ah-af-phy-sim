function [ decoder_output ] = fec_decoder( decoder_input , rate)
%Function decodes deinterleaved data input using convolution decoder K=7
%decoder_input - [bits] convolutionally coded data input [array]
%rate - [string] decoder rate 1/2, 2/3, 3/4 or 5/6 for IEEE 802.11af
%decoder_output - [B] scramled data output array

K=7; % decoder length
G0=133; % generator polynomial G0
G1=171; % generator polynomial G1
t=poly2trellis(K,[G1 G0]); % define trellis structure

data_length = length(decoder_input);
tblen_max = data_length * str2num(rate); % max traceback depth length

switch(rate)
    case '1/2' 
        rate_vector = [1;1]; %R=1/2
        tblen = 5*(K-1); % traceback depth
    case '2/3' 
        rate_vector = [1;1;1;0]; %R=2/3
        tblen = 7.5*(K-1); % traceback depth
    case '3/4' 
        rate_vector = [1;1;1;0;0;1]; %R=3/4
        tblen = 10*(K-1); % traceback depth
    case '5/6' 
        rate_vector = [1;1;1;0;0;1;1;0;0;1]; %R=5/6
        tblen = 15*(K-1); % traceback depth
    otherwise error('Coding rate incorrect!')
end

if ( tblen > tblen_max )
    tblen = tblen_max;
end

%tbl = 24;
%tbl=5*(K-1); % traceback depth
%tbl=96;
%vitdec_output = vitdec(decoder_input',t,tbl,'term','hard',rate_vector);

% if gpuDeviceCount == 0 % gpu device not available
    H = comm.ViterbiDecoder('TrellisStructure',t,...
        'InputFormat','hard',...
        'TerminationMethod','Truncated',...
        'TracebackDepth',tblen,...
        'PuncturePatternSource','Property',...
        'PuncturePattern',rate_vector);
% else
%     H = comm.gpu.ViterbiDecoder('TrellisStructure',t,...
%         'InputFormat','hard',...
%         'TerminationMethod','Truncated',...
%         'TracebackDepth',tblen,...
%         'PuncturePatternSource','Property',...
%         'PuncturePattern',rate_vector,...
%         'NumFrames',1);
% end

                        
vitdec_output = step(H, decoder_input');

decoder_output = vitdec_output';
% bits_len = length(vitdec_output);
% col=8;
% row=bits_len/col;
% decoder_output = zeros(row,col);
% 
% k=0;
% for j=1:row
%     for i=1:col
%         decoder_output(j,i) = vitdec_output(k*col+i);
%     end
%     k=k+1;
% end
% 
% decoder_output = bi2de(decoder_output,'left-msb');
% 
% decoder_output = decoder_output';

end
