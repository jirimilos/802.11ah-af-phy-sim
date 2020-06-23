function [Nsd, Nsa, Nbpsc, Ncbps, Ndbps, Kmod, pilots ] = mcss(channelWidth, param, model )

% This function, according to IEEE Std 802.11ah, associate all basic MCS parameters related to signal bandwidth

switch channelWidth
 %% channel width 1 MHz
    case 1
        Nsd = 24;
        Nsa = 32; 
        pilots = [-7, 7];
        
        switch model
            case 'SISO'
                switch (param)
                    case 0
                        Nbpsc = 1;          % coded bits per subcarrier
                        Ncbps = 24;         % coded bits per OFDM symbol
                        Ndbps = 12;         % data bits per OFDM symbol
                        Kmod  = 1;          % normalization factor
                    case 1
                        Nbpsc = 2;
                        Ncbps = 48;
                        Ndbps = 24; 
                        Kmod  = 1/sqrt(2);
                    case 2
                        Nbpsc = 2;
                        Ncbps = 48;
                        Ndbps = 36; 
                        Kmod  = 1/sqrt(2);
                    case 3 
                        Nbpsc = 4;
                        Ncbps = 96;
                        Ndbps = 48; 
                        Kmod  = 1/sqrt(10);
                    case 4
                        Nbpsc = 4;
                        Ncbps = 96;
                        Ndbps = 72; 
                        Kmod  = 1/sqrt(10);                
                    case 5
                        Nbpsc = 6;
                        Ncbps = 144;
                        Ndbps = 96; 
                        Kmod  = 1/sqrt(42);             
                    case 6
                        Nbpsc = 6;
                        Ncbps = 144;
                        Ndbps = 108; 
                        Kmod  = 1/sqrt(42); 
                    case 7
                        Nbpsc = 6;
                        Ncbps = 144;
                        Ndbps = 120; 
                        Kmod  = 1/sqrt(42); 
                    case 8
                        Nbpsc = 8;
                        Ncbps = 192;
                        Ndbps = 144; 
                        Kmod  = 1/sqrt(170);      
                    case 9
                        Nbpsc = 8;
                        Ncbps = 192;
                        Ndbps = 160; 
                        Kmod  = 1/sqrt(170); 
                    case 10
                        Nbpsc = 1;
                        Ncbps = 24;
                        Ndbps = 6; 
                        Kmod  = 1;                 
                    otherwise
                        error(['Unknown param ' num2str(param)]);
                end
            case 'MIMO'
                switch (param)
                    case 0
                        Nbpsc = 1;          % coded bits per subcarrier
                        Ncbps = 24;         % coded bits per OFDM symbol
                        Ndbps = 12;         % data bits per OFDM symbol
                        Kmod  = 1;          % normalization factor
                    case 1
                        Nbpsc = 2;
                        Ncbps = 48;
                        Ndbps = 24; 
                        Kmod  = 1/sqrt(2);
                    case 2
                        Nbpsc = 2;
                        Ncbps = 48;
                        Ndbps = 36; 
                        Kmod  = 1/sqrt(2);            
                    otherwise
                        error(['Unknown param ' num2str(param) ' for MIMO']);
                end 
        end
        
%% channel width 2 MHz
    case 2
        Nsd = 52;
        Nsa = 64;
        pilots = [-21, -7, 7, 21];
        switch model
            case 'SISO'
                switch (param)
                    case 0
                        Nbpsc = 1;          % coded bits per subcarrier
                        Ncbps = 52;         % coded bits per OFDM symbol
                        Ndbps = 26;         % data bits per OFDM symbol
                        Kmod  = 1;          % normalization factor
                    case 1
                        Nbpsc = 2;
                        Ncbps = 104;
                        Ndbps = 52; 
                        Kmod  = 1/sqrt(2);
                    case 2
                        Nbpsc = 2;
                        Ncbps = 104;
                        Ndbps = 78; 
                        Kmod  = 1/sqrt(2);
                    case 3 
                        Nbpsc = 4;
                        Ncbps = 208;
                        Ndbps = 104; 
                        Kmod  = 1/sqrt(10);
                    case 4
                        Nbpsc = 4;
                        Ncbps = 208;
                        Ndbps = 156; 
                        Kmod  = 1/sqrt(10);                
                    case 5
                        Nbpsc = 6;
                        Ncbps = 312;
                        Ndbps = 208; 
                        Kmod  = 1/sqrt(42);             
                    case 6
                        Nbpsc = 6;
                        Ncbps = 312;
                        Ndbps = 234; 
                        Kmod  = 1/sqrt(42); 
                    case 7
                        Nbpsc = 6;
                        Ncbps = 312;
                        Ndbps = 260; 
                        Kmod  = 1/sqrt(42); 
                    case 8
                        Nbpsc = 8;
                        Ncbps = 416;
                        Ndbps = 312; 
                        Kmod  = 1/sqrt(170);                 
                    otherwise
                        error(['Unknown param ' num2str(param)]);
                end
            case 'MIMO'
                switch (param)
                    case 0
                        Nbpsc = 1;       	% coded bits per subcarrier
                        Ncbps = 52;         % coded bits per OFDM symbol
                        Ndbps = 26;         % data bits per OFDM symbol
                        Kmod  = 1;          % normalization factor
                    case 1
                        Nbpsc = 2;
                        Ncbps = 104;
                        Ndbps = 52; 
                        Kmod  = 1/sqrt(2);
                    case 2
                        Nbpsc = 2;
                        Ncbps = 104;
                        Ndbps = 78; 
                        Kmod  = 1/sqrt(2);              
                    otherwise
                        error(['Unknown param ' num2str(param) ' for MIMO']);
                end 
        end
        
%% channel width 4 MHz
    case 4
        Nsd = 108;
        Nsa = 128;
        pilots = [-53, -25, -11, 11, 25, 53];
        switch model
            case 'SISO'
                switch (param)
                    case 0
                        Nbpsc = 1;           % coded bits per subcarrier
                        Ncbps = 108;         % coded bits per OFDM symbol
                        Ndbps = 54;          % data bits per OFDM symbol
                        Kmod  = 1;           % normalization factor
                    case 1
                        Nbpsc = 2;
                        Ncbps = 216;
                        Ndbps = 108; 
                        Kmod  = 1/sqrt(2);
                    case 2
                        Nbpsc = 2;
                        Ncbps = 216;
                        Ndbps = 162; 
                        Kmod  = 1/sqrt(2);
                    case 3 
                        Nbpsc = 4;
                        Ncbps = 432;
                        Ndbps = 216; 
                        Kmod  = 1/sqrt(10);
                    case 4
                        Nbpsc = 4;
                        Ncbps = 432;
                        Ndbps = 324; 
                        Kmod  = 1/sqrt(10);                
                    case 5
                        Nbpsc = 6;
                        Ncbps = 648;
                        Ndbps = 432; 
                        Kmod  = 1/sqrt(42);             
                    case 6
                        Nbpsc = 6;
                        Ncbps = 648;
                        Ndbps = 486; 
                        Kmod  = 1/sqrt(42); 
                    case 7
                        Nbpsc = 6;
                        Ncbps = 648;
                        Ndbps = 540; 
                        Kmod  = 1/sqrt(42); 
                    case 8
                        Nbpsc = 8;
                        Ncbps = 864;
                        Ndbps = 648; 
                        Kmod  = 1/sqrt(170);   
                    case 9
                        Nbpsc = 8;
                        Ncbps = 864;
                        Ndbps = 720; 
                        Kmod  = 1/sqrt(170); 
                    otherwise
                        error(['Unknown param ' num2str(param)]);
                end    
            case 'MIMO'
                switch (param)
                    case 0
                        Nbpsc = 1;           % coded bits per subcarrier
                        Ncbps = 108;         % coded bits per OFDM symbol
                        Ndbps = 54;          % data bits per OFDM symbol
                        Kmod  = 1;           % normalization factor
                    case 1
                        Nbpsc = 2;
                        Ncbps = 216;
                        Ndbps = 108; 
                        Kmod  = 1/sqrt(2);
                    case 2
                        Nbpsc = 2;
                        Ncbps = 216;
                        Ndbps = 162; 
                        Kmod  = 1/sqrt(2);            
                    otherwise
                        error(['Unknown param ' num2str(param) ' for MIMO']);
                end 
        end
        
%% channel width 8 MHz
    case 8
        Nsd = 234;
        Nsa = 256;
        pilots = [-103, -75, -39, -11, 11, 39, 75, 103];
        switch model
            case 'SISO'
                switch (param)
                    case 0
                        Nbpsc = 1;           % coded bits per subcarrier
                        Ncbps = 234;         % coded bits per OFDM symbol
                        Ndbps = 117;         % data bits per OFDM symbol
                        Kmod  = 1;           % normalization factor
                    case 1
                        Nbpsc = 2;
                        Ncbps = 468;
                        Ndbps = 234; 
                        Kmod  = 1/sqrt(2);
                    case 2
                        Nbpsc = 2;
                        Ncbps = 468;
                        Ndbps = 351; 
                        Kmod  = 1/sqrt(2);
                    case 3 
                        Nbpsc = 4;
                        Ncbps = 936;
                        Ndbps = 468; 
                        Kmod  = 1/sqrt(10);
                    case 4
                        Nbpsc = 4;
                        Ncbps = 936;
                        Ndbps = 702; 
                        Kmod  = 1/sqrt(10);                
                    case 5
                        Nbpsc = 6;
                        Ncbps = 1404;
                        Ndbps = 936; 
                        Kmod  = 1/sqrt(42);             
                    case 6
                        Nbpsc = 6;
                        Ncbps = 1404;
                        Ndbps = 1053; 
                        Kmod  = 1/sqrt(42); 
                    case 7
                        Nbpsc = 6;
                        Ncbps = 1404;
                        Ndbps = 1170; 
                        Kmod  = 1/sqrt(42); 
                    case 8
                        Nbpsc = 8;
                        Ncbps = 1872;
                        Ndbps = 1404; 
                        Kmod  = 1/sqrt(170);   
                    case 9
                        Nbpsc = 8;
                        Ncbps = 1872;
                        Ndbps = 1560; 
                        Kmod  = 1/sqrt(170); 
                    otherwise
                        error(['Unknown param ' num2str(param)]);
                end 
            case 'MIMO'
                switch (param)
                    case 0
                        Nbpsc = 1;           % coded bits per subcarrier
                        Ncbps = 234;         % coded bits per OFDM symbol
                        Ndbps = 117;         % data bits per OFDM symbol
                        Kmod  = 1;           % normalization factor
                    case 1
                        Nbpsc = 2;
                        Ncbps = 468;
                        Ndbps = 234; 
                        Kmod  = 1/sqrt(2);
                    case 2
                        Nbpsc = 2;
                        Ncbps = 468;
                        Ndbps = 351; 
                        Kmod  = 1/sqrt(2);             
                    otherwise
                        error(['Unknown param ' num2str(param) ' for MIMO']);
                end 
        end
        
%% channel width 16 MHz
    case 16
        Nsd = 468;
        Nsa = 512;
        pilots = [-231,-203,-167,-139,-117,-89,-53,-25,25,53,89,117,139,167,203,231];
        switch model
            case 'SISO'
                switch (param)
                    case 0
                        Nbpsc = 1;           % coded bits per subcarrier
                        Ncbps = 468;         % coded bits per OFDM symbol
                        Ndbps = 234;         % data bits per OFDM symbol
                        Kmod  = 1;           % normalization factor
                    case 1
                        Nbpsc = 2;
                        Ncbps = 936;
                        Ndbps = 468; 
                        Kmod  = 1/sqrt(2);
                    case 2
                        Nbpsc = 2;
                        Ncbps = 936;
                        Ndbps = 702; 
                        Kmod  = 1/sqrt(2);
                    case 3 
                        Nbpsc = 4;
                        Ncbps = 1872;
                        Ndbps = 936; 
                        Kmod  = 1/sqrt(10);
                    case 4
                        Nbpsc = 4;
                        Ncbps = 1872;
                        Ndbps = 1404; 
                        Kmod  = 1/sqrt(10);                
                    case 5
                        Nbpsc = 6;
                        Ncbps = 2808;
                        Ndbps = 1872; 
                        Kmod  = 1/sqrt(42);             
                    case 6
                        Nbpsc = 6;
                        Ncbps = 2808;
                        Ndbps = 2106; 
                        Kmod  = 1/sqrt(42); 
                    case 7
                        Nbpsc = 6;
                        Ncbps = 2808;
                        Ndbps = 2340; 
                        Kmod  = 1/sqrt(42); 
                    case 8
                        Nbpsc = 8;
                        Ncbps = 3744;
                        Ndbps = 2808; 
                        Kmod  = 1/sqrt(170);   
                    case 9
                        Nbpsc = 8;
                        Ncbps = 3744;
                        Ndbps = 3120; 
                        Kmod  = 1/sqrt(170); 
                    otherwise
                        error(['Unknown param ' num2str(param)]);
                end 
            case 'MIMO'
                switch (param)
                    case 0
                        Nbpsc = 1;           % coded bits per subcarrier
                        Ncbps = 468;         % coded bits per OFDM symbol
                        Ndbps = 234;         % data bits per OFDM symbol
                        Kmod  = 1;           % normalization factor
                    case 1
                        Nbpsc = 2;
                        Ncbps = 936;
                        Ndbps = 468; 
                        Kmod  = 1/sqrt(2);
                    case 2
                        Nbpsc = 2;
                        Ncbps = 936;
                        Ndbps = 702; 
                        Kmod  = 1/sqrt(2);           
                    otherwise
                        error(['Unknown param ' num2str(param) ' for MIMO']);
                end 
        end
        
    otherwise
        
        error(['Unknown channel width ' num2str(channelWidth)]);
end

end

