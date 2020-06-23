function [ data_Rx_equalized ] = equalizerZF( data_Rx, pilot_indexes, CH_BANDWIDTH, transmission_type )
%Zero forcing equalizer for fading channels
% data_Rx - [matrix] received OFDM data without CP
% pilot_indexes - [array] pilot indexes per frame

[N_SYM, FFT_size] = size(data_Rx);

data_Rx_equalized = zeros(N_SYM, FFT_size);
data_Rx_approx = zeros(N_SYM, FFT_size);
%H_pilots = zeros(N_SYM, size(pilot_indexes));
pilots = pilots_generator(N_SYM, CH_BANDWIDTH, transmission_type); % generates pilots structure (works only for non-HT dup mode TVHT_W!)

pilots_Rx = data_Rx(:,pilot_indexes); % pilots received

xq = 1:FFT_size; % approximation range

H_pilots = pilots_Rx./pilots; % pilot characteristic

for i=1:N_SYM
    data_Rx_approx(i,:) = interp1(pilot_indexes,H_pilots(i,:),xq,'linear','extrap');
end

%H_chan = data_Rx./data_Rx_approx;

%G_eq = 1./H_chan;

G_eq = 1./data_Rx_approx;

data_Rx_equalized = data_Rx.*G_eq;

%H_pilots = pilots_Rx./pilots; % channel characteristic

%G_pilots = 1./H_pilots; % calculated equalizer coefficients

% equalized data Rx
% for j = 1:N_SYM
%    data_Rx_equalized(j,:) = data_Rx(j,:) * G_pilots(j,1);
% end

% for j = 1:N_SYM
%     for i = 1:FFT_size
%         if (i <= 12) 
%             data_Rx_equalized(j,i) = data_Rx(j,i) * G_pilots(j,1);
%         elseif (i <= 26)
%             data_Rx_equalized(j,i) = data_Rx(j,i) * G_pilots(j,2);
%         elseif (i <= 40)
%             data_Rx_equalized(j,i) = data_Rx(j,i) * G_pilots(j,3);
%         elseif (i <= 65)
%             data_Rx_equalized(j,i) = data_Rx(j,i) * G_pilots(j,4);
%         elseif (i <= 90)
%             data_Rx_equalized(j,i) = data_Rx(j,i) * G_pilots(j,5);
%         elseif (i <= 104)
%             data_Rx_equalized(j,i) = data_Rx(j,i) * G_pilots(j,6);
%         elseif (i <= 118)
%             data_Rx_equalized(j,i) = data_Rx(j,i) * G_pilots(j,7);
%         elseif (i <= 128)
%             data_Rx_equalized(j,i) = data_Rx(j,i) * G_pilots(j,8);
%         end
%     end
% end

end

