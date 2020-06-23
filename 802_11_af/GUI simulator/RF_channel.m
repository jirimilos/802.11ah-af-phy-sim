function [ RF_signal ] = RF_channel( RF_channel_model, Tx_signal, SNR, Signal_Power, fs )
% Function process transmitted signal trough RF channel model
% RF_channel_model - AWGN/Rayleign/Ricien chosen model of RF channel
% SNR - S/N ration for AWGN channel
% fs - symbol frequency [Hz]

switch(RF_channel_model)
    case 'AWGN' % AWGN channel
        
        %SNR = 20; % 20 dB
        %powerdBW = 10 * log10(sum(abs(output_signal(:)).^2)/length(output_signal(:)));  % calculate power of output signal
        %RF_signal = awgn(Tx_signal, SNR, 'measured'); % add AWGN noise to source signal
        RF_signal = awgn(Tx_signal, SNR, Signal_Power); % add AWGN noise to source signal
    
    case 'Rician (12 paths)' % Rice fading model (RA12)
        
        K = 10; % ratio of signal power in dominant component over scattered power [dB]
        
%         ro = [7.9 10.5 11.8 14.8 8.0 ... 
%               10.7 9.2 11.7 13.0 12.5 ...
%               13.2]; % gain in each of the three paths [dB]
%         ro = 10.^(-ro/10); % convert gain to [-]
                                
        ro = [0.407163 0.303585 0.258782 0.185074 0.400967 ... 
              0.295723 0.350825 0.262909 0.225894 0.240140 ...
              0.221155 0.259730]; % gain in each of the three paths [-]        
        
        tau = [0.518650 2.751772 0.602895 3.324866 1.935570 ...
               0.429948 3.228872 0.848831 0.073883 0.924450 ...
               0.640512 1.368671]; % vector of delay [us]

        theta = (pi/180).*([336.0 127.0 215.3 330.9 8.8 ...
                            339.7 174.9 36.0 122.0 210.0 ...
                            191.0 22.6]); % vector of phase shift [rad]
                        
        % Gain calculation for direct path according to K and adding values
        % of direct path to vectors ro, tau and theta
        ro0 = sqrt(K*sum(ro.^2)); % define power on the direct path
        ro = [ro0 ro]; % adding power of the direct path
        tau = [0 tau].*1E-6; % adding delay of the direct path [s]
        theta = [0 theta]; % adding phase shift of the direct path
        
        %Ts = (1/fs); % one ofdm symbol period
        Ts = 3.2e-6;
        
        n_paths = length(ro); % number of paths
        signal_length = length(Tx_signal); % signal length
        signal_path = zeros(1,signal_length); % signal path
        
        RF_signal = zeros(1,signal_length); % signal length with delay
        
        % main loop for delay calculations
        for path = 1:n_paths 
            signal_path(1,:) = Tx_signal;
            signal_path(1,:) = ro(1,path)*exp(-1i*theta(path))*signal_path(1,:);
            delay = round(tau(1,path)/Ts);
            signal_path(1,:) = [zeros(1,delay) signal_path(1,1:signal_length-delay)];
            RF_signal = RF_signal + signal_path(1,:);
        end
                
        % normalization
        norm = sqrt(sum(ro.^2));
        RF_signal = RF_signal./norm;
        
%         pz_gain = 10.^(-ro/20);
%         shifting = fix(tau./Ts);
%         Max_shifting = max(shifting);
%         RF_signal = zeros(1,signal_length+Max_shifting);
%         
%         % Modeling of the Ricean channel
%         for path = 1:n_paths
%             Signal_Path = zeros(1,signal_length+Max_shifting);
%             Signal_Path(shifting(path)+1:shifting(path)+signal_length) = Tx_signal;
%             gain = pz_gain(path)*exp(-1i*theta(path));
%             Signal_Path = Signal_Path*gain;
%             RF_signal = RF_signal + Signal_Path;
%         end
%         
%         RF_signal = RF_signal(1:signal_length);
%         norm = sqrt(sum(ro.^2));
%         RF_signal = RF_signal/norm;
        
%         coeff = 1.9; 
%         RF_signal = awgn(RF_signal, SNR, coeff*Signal_Power); % add AWGN noise to source signal
        RF_signal = awgn(RF_signal, SNR, Signal_Power); % add AWGN noise to source signal

    case 'Rayleigh (12 paths)' % Rayleigh model (no direct path)
        
%         ro = [7.9 10.5 11.8 14.8 8.0 ... 
%               10.7 9.2 11.7 13.0 12.5 ...
%               13.2 11.8]; % gain in each of the three paths [dB]
%         ro = 10.^(-ro/10); % convert gain to [-]
%         tau = [0.5 2.75 0.6 3.3 1.95 ...
%                0.45 3.25 0.85 0.05 0.9 ...
%                0.65 1.35]; % vector of delay [us]

        ro = [0.407163 0.303585 0.258782 0.185074 0.400967 ... 
              0.295723 0.350825 0.262909 0.225894 0.240140 ...
              0.221155 0.259730]; % gain in each of the three paths [-]        
        
        tau = [0.518650 2.751772 0.602895 3.324866 1.935570 ...
               0.429948 3.228872 0.848831 0.073883 0.924450 ...
               0.640512 1.368671]; % vector of delay [us]

        theta = (pi/180).*([336.0 127.0 215.3 330.9 8.8 ...
                            339.7 174.9 36.0 122.0 210.0 ...
                            191.0 22.6]); % vector of phase shift [rad]
        
        %Ts = (1/fs); % one symbol sampling period
        Ts = 2.1e-6;
        tau = tau.*1E-6; % delay conversion to seconds
        %rxo=ro; % power on the each path
        %ro = 10.^(-rxo/20); % powers on the each path [-]
        n_paths = length(ro); % number of paths
        signal_length = length(Tx_signal); % signal length
        signal_path = zeros(1,signal_length); % signal path
        RF_signal = zeros(1,signal_length); % signal length with delay
        
       % main loop for delay calculations
        for path = 1:n_paths 
            signal_path(1,:) = Tx_signal;
            signal_path(1,:) = ro(1,path)*exp(-1i*theta(path))*signal_path(1,:);
            delay = round(tau(1,path)/Ts);
            signal_path(1,:) = [zeros(1,delay) signal_path(1,1:signal_length-delay)];
            RF_signal = RF_signal + signal_path(1,:);
         end
        
        % normalization
        norm = sqrt(sum(ro.^2));
        RF_signal = RF_signal./norm;
        
        %coeff = 1; 
        %RF_signal = awgn(RF_signal, SNR, coeff*Signal_Power); % add AWGN noise to source signal
        RF_signal = awgn(RF_signal, SNR, Signal_Power); % add AWGN noise to source signal        
   
    case 'PI' % Pedestrian Indoor model for 3km/h receiver speed
        
        K = 10; % Rician factor [dB]
        
        ro = [0.0 -6.4 -10.4 -13.0 -13.3 ...
            -13.7 -16.2 -15.2 -14.9 -16.2 ...
            -11.1 -11.2]; % gain in each of the 12 paths [dB]
        
        tau = [0.0 0.1 0.2 0.4 0.6 ...
                0.8 1.0 1.6 8.1 8.8 ...
                9.0 9.2]; % delay in each of the 12 paths [us] 
        
        T_OFDM = 1/(1*fs); % OFDM symbol period [s]
        T_GI = 4.5e-6; % GI duration [s]
        Ts = T_OFDM + T_GI; % one symbol period [s]
        %Ts = 100e-6;
        ts = 1e-10;
        %fc=626;
        fc = 868;   % carrier frequency [MHz]
        N = 15; % oversampling factor
        Fs = N *( fc * (10 ^ 6) ); % sampling frequency [Hz]
        %Fs = 1/ts;
        vc = 3e8;   % speed of light  [m/s]
        v = 3; % speed of pedestrian [km/h]
        fd = v*fc/(3.6*300); % maximum Doppler shift [Hz]
        %tau = tau.*1e-6; % conversion to seconds      
        %fd_ratio = 0.5; % Doppler frequency ratio
        fd_ratio = fd * 0.5; % Doppler frequency ratio
        STD_norm = 0.08; % normalized deviation STD
        
        ro_gain = 10.^(ro/20); % convert gain from dB [-]
        n_paths = length(tau); % number of paths
        signal_length = length(Tx_signal); % length of signal

        RF_signal = zeros(1,length(Tx_signal));
        signal_path = zeros(n_paths,signal_length);

%         % main loop for delay calculations         
%         for path = 1:n_paths
%          
%             if ro_gain(1,path)==1              % rice - gauss
%                
%                 ch = ricianchan(ts,fd,K);
%                 ch.DirectPathDopplerShift = fd_ratio*fd;
%                 ch.NormalizePathGains = 0;                          
%                 ch.ResetBeforeFiltering = 1;
%                 ch.DopplerSpectrum = doppler.gaussian(STD_norm);
%    
%                 signal_path(path,:) = filter(ch,Tx_signal);
%                 signal_path(path,:) = ro_gain(1,path)*signal_path(path,:);
%                 del = round(tau(1,path)*1E-6/Ts);
%                 signal_path(path,:) = [zeros(1,del) signal_path(path,1:signal_length-del)];
%                 RF_signal = RF_signal + signal_path(path,:);
%                        
%                           
%             else                               % rayleigh - gauss
%                  
%                 ch = rayleighchan(ts,fd);  
%                 ch.NormalizePathGains = 0;                          
%                 ch.ResetBeforeFiltering = 1;
%                 ch.DopplerSpectrum = doppler.gaussian(STD_norm);
%                     
%                 signal_path(path,:) = filter(ch,Tx_signal);
%                 signal_path(path,:) = ro_gain(1,path)*signal_path(path,:);
%                 del = round(tau(1,path)*1E-6/Ts);
%                 signal_path(path,:) = [zeros(1,del) signal_path(path,1:signal_length-del)];
%                 RF_signal = RF_signal + signal_path(path,:);
%                 
%             end
%              
%          end
%         % normalization
%         norm = sqrt(sum(ro_gain.^2));
%         RF_signal = RF_signal./norm;

%         % main loop for delay calculations         
%         for path = 1:n_paths
%          
%             if ro_gain(1,path)==1              % rice - gauss
%                
%                 ch = ricianchan(ts,fd,K);
%                 ch.DirectPathDopplerShift = fd_ratio*fd;
%                 ch.NormalizePathGains = 0;                          
%                 ch.ResetBeforeFiltering = 1;
%                 ch.DopplerSpectrum = doppler.gaussian(STD_norm);
%    
%                 signal_path(1,:) = filter(ch,Tx_signal);
%                 signal_path(1,:) = ro_gain(1,path)*signal_path(1,:);
%                 del = round(tau(1,path)*1E-6/Ts);
%                 signal_path(1,:) = [zeros(1,del) signal_path(1,1:signal_length-del)];
%                 RF_signal = RF_signal + signal_path(1,:);
%                        
%                           
%             else                               % rayleigh - gauss
%                  
%                 ch = rayleighchan(ts,fd);  
%                 ch.NormalizePathGains = 0;                          
%                 ch.ResetBeforeFiltering = 1;
%                 ch.DopplerSpectrum = doppler.gaussian(STD_norm);
%                     
%                 signal_path(1,:) = filter(ch,Tx_signal);
%                 signal_path(1,:) = ro_gain(1,path)*signal_path(1,:);
%                 del = round(tau(1,path)*1E-6/Ts);
%                 signal_path(1,:) = [zeros(1,del) signal_path(1,1:signal_length-del)];
%                 RF_signal = RF_signal + signal_path(1,:);
%                 
%             end
%              
%         end
%         

        rayChan = comm.RayleighChannel(...
            'SampleRate',Fs,...
            'NormalizePathGains',1,...
            'RandomStream','mt19937ar with seed', ...
            'Seed',73, ...
            'MaximumDopplerShift',fd,...
            'DopplerSpectrum',doppler('Gaussian',STD_norm),...
            'PathGainsOutputPort',true);
        
        ricianChan = comm.RicianChannel(...
            'SampleRate',Fs,...
            'NormalizePathGains',1,...
            'RandomStream','mt19937ar with seed', ...
            'Seed',73, ...
            'KFactor',K,...
            'DirectPathDopplerShift',(fd_ratio),...
            'MaximumDopplerShift',fd,...
            'DopplerSpectrum',doppler('Gaussian',STD_norm),...
            'PathGainsOutputPort',true);

        RF_signal = zeros(length(Tx_signal),1);
        signal_path = zeros(signal_length,n_paths);
        
        reset(rayChan);
        reset(ricianChan);
        
        %main loop for delay calculations         
        for path = 1:n_paths
            
            if ( ro_gain(1,path) == 1 )% Rician - Gauss

                signal_path(:,path) = step(ricianChan, Tx_signal');
                signal_path(:,path) = ro_gain(1,path)*signal_path(:,path);
                del = round( ( tau(1,path) * 1E-6 ) / Ts);
                signal_path(:,path) = [zeros(del,1) ; signal_path(1:signal_length-del,path)];
                RF_signal = RF_signal + signal_path(:,path);  
                
            else % Rayleigh - Gauss

                signal_path(:,path) = step(rayChan, Tx_signal');
                signal_path(:,path) = ro_gain(1,path)*signal_path(:,path);
                del = round( ( tau(1,path) * 1E-6 ) / Ts);
                signal_path(:,path) = [zeros(del,1) ; signal_path(1:signal_length-del,path)];
                RF_signal = RF_signal + signal_path(:,path);            
            end
         
        end
        
        RF_signal_tmp = RF_signal;
        RF_signal = zeros(1,signal_length);
        
        for i=1:signal_length
            RF_signal(i) = RF_signal_tmp(i,1);
        end
        
        %normalization
        norm = sqrt(sum(ro_gain.^2));
        RF_signal = RF_signal./norm;
        
        RF_signal = awgn(RF_signal, SNR, Signal_Power); % add AWGN noise to source signal

    otherwise % No fading
        
        %RF_signal = output_signal;
        error('Choosen RF channel model does not exist!')
end


end

