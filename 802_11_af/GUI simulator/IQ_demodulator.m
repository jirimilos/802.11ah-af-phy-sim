function [ I, Q ] = IQ_demodulator(RF_signal, carrier_I, carrier_Q)
% Function process IQ demodulation for received RF signal
% RF_signal - Received from communication channel RF signal
% carrier_I - Real part of the frequency carrier
% carrier_Q - Imaginary part of the frequency carrier

I = RF_signal.*carrier_I;
Q = RF_signal.*carrier_Q;

end

