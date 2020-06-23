function [ ekvalizChan ] = equalizer( fftData, pilots, refPilots, Nsa)
% This function provides Zero-Forcing (ZF) equalization
% In the case of AWGN channel (no ISI), this function is not used  

        % Size of parameters
        [nRow, nCol] = size(refPilots);
        chanPilots = zeros(nRow, nCol);
        
        % Selecting pilots after going through a channel 
        for i = 1:length(pilots)
            chanPilots(i,:) = fftData(pilots(i),:);
        end
               
        xq = 1:Nsa; % for QPSK 1:64
        
        % Least Square estimation (LS)
        
        hat_pilots_LS = zeros(nRow,nCol);
        H_LS_interp = zeros(Nsa,nCol);

        for ii = 1:size(refPilots,2)
            hat_pilots_LS(:,ii) = chanPilots(:,ii)./refPilots(:,ii);
            H_LS_interp(:,ii) = interp1(pilots,hat_pilots_LS(:,ii),xq,'linear','extrap').';
        end
        
        bb = sum(H_LS_interp,2)./nCol;
           
        for subcarrier = 1:size(H_LS_interp,1)
            inv_tmp(subcarrier,1) = pinv(bb(subcarrier,1));
        end
        
        cc = repmat(inv_tmp,1,nCol);
        ekvalizChan = cc.*fftData; 

end

