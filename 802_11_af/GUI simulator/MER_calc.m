function [ MER_dB ] = MER_calc( Tx_mod, Rx_mod )
% Function calculates Modulation Error Ratio
% MER_dB - MER calculated in dB
% Tx_mod - Modulated symbols in transmitter (referenced)
% Rx_mod - Modulated symbols in receiver (measured)

N = length(Tx_mod); % number of symbols 

I_ref = real(Tx_mod); % reference vector of I components
Q_ref = imag(Tx_mod); % reference vector of Q components

I_meas = real(Rx_mod); % measured vector of I components
Q_meas = imag(Rx_mod); % measured vector of Q components

%IQ_meas = 0; % sum of measured I^2+Q^2 components
IQ_ref = 0; % sum of reference I^2+Q^2 components
IQ_diff = 0; % sum of differences between measured and referenced I and Q components

for i=1:N
    
    %IQ_meas_act = I_meas(i) ^ 2 + Q_meas(i) ^ 2;
    %IQ_meas = IQ_meas + IQ_meas_act;
    IQ_ref_act = I_ref(i) ^ 2 + Q_ref(i) ^ 2;
    IQ_ref = IQ_ref + IQ_ref_act;
    
    IQ_diff_act = ( I_ref(i) - I_meas(i) ) ^ 2 + ( Q_ref(i) - Q_meas(i) ) ^ 2;
    IQ_diff = IQ_diff + IQ_diff_act;
    
end

%MER_dB = 10 * log10 ( IQ_meas / IQ_diff );
MER_dB = 10 * log10 ( IQ_ref / IQ_diff );

end

