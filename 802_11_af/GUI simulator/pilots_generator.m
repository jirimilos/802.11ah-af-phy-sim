function [ pilots ] = pilots_generator( N_SYM, CH_BANDWIDTH,FORMAT )
% Function generates pilots (necessary for receiver)
% implemented for non-HT duplicate and VHT modes
% TVHT_W

% N_SYM - [int] number of data OFDM symbols 
% FORMAT - [char] transmission type DUP Non-HT or VHT

switch(CH_BANDWIDTH)
    case 'TVHT_W'
        IFFT_size = 128;
        switch(FORMAT)
            case 'non-HT'
                subcarrier_non_HT_indexes = [-26:-1,1:26]; % subcarrier indexes per non-HT mode
                subcarriers_DC = zeros(1,(10+1)); % subcarriers (1 DC + 10 NULLs)

                subcarriers_non_HT = [ % subcarriers for non-HT mode
                    subcarrier_non_HT_indexes(find(subcarrier_non_HT_indexes == -26):find(subcarrier_non_HT_indexes == -1)) 0 ... 
                    subcarrier_non_HT_indexes(find(subcarrier_non_HT_indexes == 1):find(subcarrier_non_HT_indexes == 26))];

                subcarriers_non_HT_dup = [ % subcarriers for non-HT duplicate mode
                    subcarriers_non_HT subcarriers_DC subcarriers_non_HT];

                % NULLs addition to fill IFFT block 
                fft_nulls =IFFT_size - length(subcarriers_non_HT_dup); % number of NULLs to fill IFFT block
                a_nulls = ceil(fft_nulls/2); % NULLs on the begin
                b_nulls = fft_nulls-a_nulls; % NULLs on the end
                subcarriers_non_HT_dup = [zeros(1,a_nulls) subcarriers_non_HT_dup zeros(1,b_nulls)]; % IFFT block mapping

                % pilots mapping

                pilot1_indexes = find( subcarriers_non_HT_dup == -21 );
                pilot2_indexes = find( subcarriers_non_HT_dup == -7 );
                pilot3_indexes = find( subcarriers_non_HT_dup == 7 );
                pilot4_indexes = find( subcarriers_non_HT_dup == 21 );

                pilot_indexes = [pilot1_indexes ... % subcarriers with index -21
                                pilot2_indexes ... % subcarriers with index -7
                                pilot3_indexes ... % subcarriers with index +7
                                pilot4_indexes];% subcarriers with index +21
                            
                        
            case 'VHT'

                subcarrier_VHT_indexes = [-58:-2,2:58]; % subcarrier indexes per VHT mode
                subcarriers_DC = zeros(1,(2+1)); % subcarriers (1 DC + 2 NULLs)

                subcarriers_VHT = [ % subcarriers for VHT mode include DC
                    subcarrier_VHT_indexes(find(subcarrier_VHT_indexes == -58):find(subcarrier_VHT_indexes == -2)) subcarriers_DC ... 
                    subcarrier_VHT_indexes(find(subcarrier_VHT_indexes == 2):find(subcarrier_VHT_indexes == 58))];

                % NULLs addition to fill IFFT block 
                fft_nulls =IFFT_size - length(subcarriers_VHT); % number of NULLs to fill IFFT block
                a_nulls = ceil(fft_nulls/2); % NULLs on the begin
                b_nulls = fft_nulls-a_nulls; % NULLs on the end
                subcarriers_VHT = [zeros(1,a_nulls) subcarriers_VHT zeros(1,b_nulls)]; % IFFT block mapping

                % pilots mapping

                pilot1_indexes = find( subcarriers_VHT == -53 );
                pilot2_indexes = find( subcarriers_VHT == -25 );
                pilot3_indexes = find( subcarriers_VHT == -11 );
                pilot4_indexes = find( subcarriers_VHT == 11 );
                pilot5_indexes = find( subcarriers_VHT == 25 );
                pilot6_indexes = find( subcarriers_VHT == 53 );

                pilot_indexes = [pilot1_indexes ... % subcarriers with index -53
                                pilot2_indexes ... % subcarriers with index -25
                                pilot3_indexes ... % subcarriers with index -11
                                pilot4_indexes ...% subcarriers with index +11
                                pilot5_indexes ... % subcarriers with index +25
                                pilot6_indexes]; % subcarriers with index +53
                N_SP = 6;
        end    
        
    case 'TVHT_2W' 
        IFFT_size = 512;
    case 'TVHT_4W'
     IFFT_size = 1024;
    case 'TVHT_W+W'
        IFFT_size = 128;
    case 'TVHT_2W+2W' 
        IFFT_size = 512;
    otherwise error('Incorrect channel bandwidth used! Could be TVHT_(W/2W/4W/W+W/2W+2W)')
end
% 
% IFFT_size = 128;
% subcarrier_non_HT_indexes = [-26:-1,1:26]; % subcarrier indexes per non-HT mode
% subcarriers_DC = zeros(1,(10+1)); % subcarriers (1 DC + 10 NULLs)
% 
% subcarriers_non_HT = [ % subcarriers for non-HT mode
%     subcarrier_non_HT_indexes(find(subcarrier_non_HT_indexes == -26):find(subcarrier_non_HT_indexes == -1)) 0 ... 
%     subcarrier_non_HT_indexes(find(subcarrier_non_HT_indexes == 1):find(subcarrier_non_HT_indexes == 26))];
% 
% subcarriers_non_HT_dup = [ % subcarriers for non-HT duplicate mode
%     subcarriers_non_HT subcarriers_DC subcarriers_non_HT];
% 
% % NULLs addition to fill IFFT block 
% fft_nulls =IFFT_size - length(subcarriers_non_HT_dup); % number of NULLs to fill IFFT block
% a_nulls = ceil(fft_nulls/2); % NULLs on the begin
% b_nulls = fft_nulls-a_nulls; % NULLs on the end
% subcarriers_non_HT_dup = [zeros(1,a_nulls) subcarriers_non_HT_dup zeros(1,b_nulls)]; % IFFT block mapping
% 
% % pilots mapping
% 
% pilot1_indexes = find( subcarriers_non_HT_dup == -21 );
% pilot2_indexes = find( subcarriers_non_HT_dup == -7 );
% pilot3_indexes = find( subcarriers_non_HT_dup == 7 );
% pilot4_indexes = find( subcarriers_non_HT_dup == 21 );
% 
% pilot_indexes = [pilot1_indexes ... % subcarriers with index -21
%                 pilot2_indexes ... % subcarriers with index -7
%                 pilot3_indexes ... % subcarriers with index +7
%                 pilot4_indexes];% subcarriers with index +21
% 

% %% pilots initialization
% % pn sequence of pilots polarity for 127 OFDM symbols
% p127 = [1,1,1,1, -1,-1,-1,1, -1,-1,-1,-1, 1,1,-1,1, -1,-1,1,1, -1,1,1,-1, 1,1,1,1, 1,1,-1,1, 1,1,-1,1, 1,-1,-1,1, 1,1,-1,1, -1,-1,-1,1, -1,1,-1,-1, 1,-1,-1,1, 1,1,1,1, -1,-1,1,1, -1,-1,1,-1, 1,-1,1,1, -1,-1,-1,1, 1,-1,-1,-1, -1,1,-1,-1, 1,-1,1,1, 1,1,-1,1, -1,1,-1,1, -1,-1,-1,-1, -1,1,-1,1, 1,-1,1,-1, 1,1,1,-1, -1,1,-1,-1, -1,1,1,1, -1,-1,-1,-1, -1,-1,-1];
% % pilots sequence of pilots -21, -7, +7, +21
% P = [1 1 1 -1];
% p127_coeff = ceil(N_SYM/127);
% if ( p127_coeff == 1)
%     pn = p127;
% else % extend pn sequency according to number of OFDM symbols
%     pn = repmat(p127,1,p127_coeff);
% end
% 
% % copy subcarriers for non-HT duplicate by N_SYM times 
% subcarriers = repmat(subcarriers_non_HT_dup,[N_SYM 1]);
% 
% for n = 1:N_SYM
%     
%     subcarriers(n,pilot1_indexes) = P(1) * pn(n);
%     subcarriers(n,pilot2_indexes) = P(2) * pn(n);
%     subcarriers(n,pilot3_indexes) = P(3) * pn(n);
%     subcarriers(n,pilot4_indexes) = P(4) * pn(n);
%     
% end
% 
% pilots = subcarriers(:,pilot_indexes);


if( strcmp(FORMAT,'non-HT') )

    %% pilots initialization
    % pn sequence of pilots polarity for 127 OFDM symbols
    p127 = [1,1,1,1, -1,-1,-1,1, -1,-1,-1,-1, 1,1,-1,1, -1,-1,1,1, -1,1,1,-1, 1,1,1,1, 1,1,-1,1, 1,1,-1,1, 1,-1,-1,1, 1,1,-1,1, -1,-1,-1,1, -1,1,-1,-1, 1,-1,-1,1, 1,1,1,1, -1,-1,1,1, -1,-1,1,-1, 1,-1,1,1, -1,-1,-1,1, 1,-1,-1,-1, -1,1,-1,-1, 1,-1,1,1, 1,1,-1,1, -1,1,-1,1, -1,-1,-1,-1, -1,1,-1,1, 1,-1,1,-1, 1,1,1,-1, -1,1,-1,-1, -1,1,1,1, -1,-1,-1,-1, -1,-1,-1];
    % pilots sequence of pilots -21, -7, +7, +21
    P = [1 1 1 -1];
    p127_coeff = ceil(N_SYM/127);
    if ( p127_coeff == 1)
        pn = p127;
    else % extend pn sequency according to number of OFDM symbols
        pn = repmat(p127,1,p127_coeff);
    end
    
    % copy subcarriers for non-HT duplicate by N_SYM times 
    subcarriers = repmat(subcarriers_non_HT_dup,[N_SYM 1]);

    for n = 1:N_SYM

        subcarriers(n,pilot1_indexes) = P(1) * pn(n);
        subcarriers(n,pilot2_indexes) = P(2) * pn(n);
        subcarriers(n,pilot3_indexes) = P(3) * pn(n);
        subcarriers(n,pilot4_indexes) = P(4) * pn(n);

    end

    pilots = subcarriers(:,pilot_indexes);
    
elseif( strcmp(FORMAT,'VHT') )
    
    %% pilots initialization
    % pilots sequence of pilots -53 -25 -11 +11 +25 +53 for 1 spatial stream
    P = [1, 1, 1, -1, -1, 1];
    
    % pilots values for N_SYM symbols  
    pilots = zeros(N_SYM, N_SP);
    
    for i = 1:N_SYM
        
        %calculate P sequence indexes for corresponding N_SYM symbol
        psi_1 = mod(i, N_SP);
        psi_2 = mod(i+1, N_SP);
        psi_3 = mod(i+2, N_SP);
        psi_4 = mod(i+3, N_SP);
        psi_5 = mod(i+4, N_SP);
        psi_6 = mod(i+5, N_SP);
        
        %pn sequence for corresponding N_SYM symbol
        pn = [P(psi_1+1) P(psi_2+1) P(psi_3+1) P(psi_4+1) P(psi_5+1) P(psi_6+1)];
        
        %fill pilots with pn sequence
        pilots(i,:) = pn;
        
    end
    
    % copy subcarriers for VHT by N_SYM times 
    subcarriers = repmat(subcarriers_VHT,[N_SYM 1]);
    
    for n = 1:N_SYM

        subcarriers(n,pilot1_indexes) =  pilots(n,1);
        subcarriers(n,pilot2_indexes) = pilots(n,2);
        subcarriers(n,pilot3_indexes) = pilots(n,3);
        subcarriers(n,pilot4_indexes) = pilots(n,4);
        subcarriers(n,pilot5_indexes) = pilots(n,5);
        subcarriers(n,pilot6_indexes) = pilots(n,6);
        
    end
    
end

end

