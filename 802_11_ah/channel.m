function [ signalTXChan ] = channel( signalTX, model, channelType )
  
switch model
    case 'SISO' % when SISO transmission mode is selected
         
        switch channelType
            case 'Rician'
                % Definition the parameters of the RC20 channel model
                k = 10; % K factor
                
                % path gain
                ro = [0.057662 0.176809 0.407163 0.303585 0.258782 ...
                      0.061831 0.150340 0.051534 0.185074 0.400967 ...
                      0.295723 0.350825 0.262909 0.225894 0.170996 ...
                      0.149723 0.240140 0.116587 0.221155 0.259730]; 
                  
                % delay [us]
                tau = [1.003019 5.422091 0.518650 2.751772 0.602895 ...
                       1.016585 0.143556 0.153832 3.324866 1.935570 ...
                       0.429948 3.228872 0.848831 0.073883 0.203952 ...
                       0.194207 0.924450 1.381320 0.640512 1.368671]; 
                   
                % phase [rad]
                theta = [4.855121 3.419109 5.864470 2.215894 3.758058 ...
                         5.430202 3.952093 1.093586 5.775198 0.154459 ...
                         5.928383 3.053023 0.628578 2.128544 1.099463 ...
                         3.462951 3.664773 2.833799 3.334290 0.393889]; 
                
                % calculation the gain of direct path according to the defined parameters of the channel model 
                ro0 = sqrt(k*sum(ro.^2));
                ro = [ro0 ro];
                tau = [0 tau].*1E-6; % recalculate to seconds
                theta = [0 theta];
                
                % the number of paths and length of the signal
                nPaths = length(ro); 
                signalLength = length(signalTX);
                signalPath = zeros(1,signalLength);
                
                % the length of the signal with delay
                signalTXChan = zeros(1,signalLength);

                % main loop for calculation of all delays
                for path = 1:nPaths 
                    signalPath(1,:) = signalTX;
                    h = ro(1,path)*exp(-1i*theta(path));
                    signalPath(1,:) = h*signalPath(1,:);
                    delay = round(tau(1,path)/3.2e-6);
                    signalPath(1,:) = [zeros(1,delay) signalPath(1,1:signalLength-delay)];
                    signalTXChan = signalTXChan + signalPath(1,:);
                end
                
                % normalization
                norm = sqrt(sum(ro.^2));
                signalTXChan = signalTXChan./norm;
                
            case 'Rayleigh'
                % path gain
                ro = [0.057662 0.176809 0.407163 0.303585 0.258782 ...
                      0.061831 0.150340 0.051534 0.185074 0.400967 ...
                      0.295723 0.350825 0.262909 0.225894 0.170996 ...
                      0.149723 0.240140 0.116587 0.221155 0.259730]; 
                  
                % delay [us]
                tau = [1.003019 5.422091 0.518650 2.751772 0.602895 ...
                       1.016585 0.143556 0.153832 3.324866 1.935570 ...
                       0.429948 3.228872 0.848831 0.073883 0.203952 ...
                       0.194207 0.924450 1.381320 0.640512 1.368671]; 
                   
                % phase [rad]
                theta = [4.855121 3.419109 5.864470 2.215894 3.758058 ...
                         5.430202 3.952093 1.093586 5.775198 0.154459 ...
                         5.928383 3.053023 0.628578 2.128544 1.099463 ...
                         3.462951 3.664773 2.833799 3.334290 0.393889]; 
                     
                % recalculate to seconds
                tau = tau.*1E-6;

                % the number of paths and length of the signal
                nPaths = length(ro); 
                signalLength = length(signalTX);
                signalPath = zeros(1,signalLength);

                % the length of the signal with delay
                signalTXChan = zeros(1,signalLength);

                % main loop for calculation of all delays
                for path = 1:nPaths 
                    signalPath(1,:) = signalTX;
                    signalPath(1,:) = ro(1,path)*exp(-1i*theta(path))*signalPath(1,:);
                    delay = round(tau(1,path)/2.1e-6);
                    signalPath(1,:) = [zeros(1,delay) signalPath(1,1:signalLength-delay)];
                    signalTXChan = signalTXChan + signalPath(1,:);
                 end

                % normalization
                norm = sqrt(sum(ro.^2));
                signalTXChan = signalTXChan./norm;
        end

    case 'MIMO' % in the case, when MIMO transmission mode is selected       
                % the process (defintion of channel model) is the same as in the case of SISO
                % but we have two paths... (a 2 x 2 MIMO scheme)
                % currently we have implemented only the RC20 channel model 
        
        switch channelType
            case 'Rician'
                k = 10; % K factor
                % path gain
                ro = [0.057662 0.176809 0.407163 0.303585 0.258782 ...
                      0.061831 0.150340 0.051534 0.185074 0.400967 ...
                      0.295723 0.350825 0.262909 0.225894 0.170996 ...
                      0.149723 0.240140 0.116587 0.221155 0.259730]; 
                % delay [us]
                tau = [1.003019 5.422091 0.518650 2.751772 0.602895 ...
                       1.016585 0.143556 0.153832 3.324866 1.935570 ...
                       0.429948 3.228872 0.848831 0.073883 0.203952 ...
                       0.194207 0.924450 1.381320 0.640512 1.368671]; 
                % phase [rad]
                theta = [4.855121 3.419109 5.864470 2.215894 3.758058 ...
                         5.430202 3.952093 1.093586 5.775198 0.154459 ...
                         5.928383 3.053023 0.628578 2.128544 1.099463 ...
                         3.462951 3.664773 2.833799 3.334290 0.393889]; 

                % calculation the gain of direct path according to K and adding direct path values to the vector ro,tau,theta 
                ro0 = sqrt(k*sum(ro.^2));
                ro = [ro0 ro];
                tau = [0 tau].*1E-6;
                theta = [0 theta];
                
                % number of paths and length of signal
                nPaths = length(ro); 
                signalLength = length(signalTX);
                signalPath = zeros(1,signalLength);
                
                % length of signal with delay
                signalTXChan11 = zeros(1,signalLength);
                signalTXChan12 = zeros(1,signalLength);

                % main loop for calculation of all delays
                for path = 1:nPaths 
                    signalPath(1,:,1) = signalTX(1,:,2);
                    h1 = ro(1,path)*exp(-1i*theta(path));
                    signalPath(1,:) = h1*signalPath(1,:);
                    delay = round(tau(1,path)/3.2e-6);
                    signalPath(1,:) = [zeros(1,delay) signalPath(1,1:signalLength-delay)];
                    signalTXChan12 = signalTXChan12 + signalPath(1,:);
                end
                
                for path = 1:nPaths 
                    signalPath(1,:,1) = signalTX(1,:,1);
                    h0 = ro(1,path)*exp(-1i*theta(path));
                    signalPath(1,:) = h0*signalPath(1,:);
                    delay = round(tau(1,path)/3.2e-6);
                    signalPath(1,:) = [zeros(1,delay) signalPath(1,1:signalLength-delay)];
                    signalTXChan11 = signalTXChan11 + signalPath(1,:);
                end
                
                % normalization
                norm = sqrt(sum(ro.^2));
                signalTXChan12 = signalTXChan12./norm;
                signalTXChan11 = signalTXChan11./norm;                
                signalTXChan(:,:,1) = signalTXChan12 + signalTXChan11;
               
                % path gain
                ro = [0.057662 0.176809 0.407163 0.303585 0.258782 ...
                      0.061831 0.150340 0.051534 0.185074 0.400967 ...
                      0.295723 0.350825 0.262909 0.225894 0.170996 ...
                      0.149723 0.240140 0.116587 0.221155 0.259730]; 
                % delay [us]
                tau = [1.003019 5.422091 0.518650 2.751772 0.602895 ...
                       1.016585 0.143556 0.153832 3.324866 1.935570 ...
                       0.429948 3.228872 0.848831 0.073883 0.203952 ...
                       0.194207 0.924450 1.381320 0.640512 1.368671]; 
                % phase [rad]
                theta = [4.855121 3.419109 5.864470 2.215894 3.758058 ...
                         5.430202 3.952093 1.093586 5.775198 0.154459 ...
                         5.928383 3.053023 0.628578 2.128544 1.099463 ...
                         3.462951 3.664773 2.833799 3.334290 0.393889]; 
                     
                % calculation the gain of direct path according to K and adding direct path values to the vector ro,tau,theta
                ro0 = sqrt(k*sum(ro.^2));
                ro = [ro0 ro];
                tau = [0 tau].*1E-6; 
                theta = [0 theta];
                
                % number of paths and length of signal
                nPaths = length(ro); 
                signalLength = length(signalTX);
                signalPath = zeros(1,signalLength);
                
                % length of signal with delay
                signalTXChan21 = zeros(1,signalLength);
                signalTXChan22 = zeros(1,signalLength);

                % main loop for calculation of all delays
                for path = 1:nPaths 
                    signalPath(1,:,1) = signalTX(1,:,1);
                    h2 = ro(1,path)*exp(-1i*theta(path));
                    signalPath(1,:) = h2*signalPath(1,:);
                    delay = round(tau(1,path)/3.2e-6);
                    signalPath(1,:) = [zeros(1,delay) signalPath(1,1:signalLength-delay)];
                    signalTXChan21 = signalTXChan21 + signalPath(1,:);
                end
                
                for path = 1:nPaths 
                    signalPath(1,:,1) = signalTX(1,:,2);
                    h3 = ro(1,path)*exp(-1i*theta(path));
                    signalPath(1,:) = h3*signalPath(1,:);
                    delay = round(tau(1,path)/3.2e-6);
                    signalPath(1,:) = [zeros(1,delay) signalPath(1,1:signalLength-delay)];
                    signalTXChan22 = signalTXChan22 + signalPath(1,:);
                end
                
                % normalization
                norm = sqrt(sum(ro.^2));
                signalTXChan21 = signalTXChan21./norm;
                signalTXChan22 = signalTXChan22./norm;                
                signalTXChan(:,:,2) = signalTXChan21 + signalTXChan22;
        end  
end

end

