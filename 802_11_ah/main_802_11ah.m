if ~(exist('usedGUI','var'))
clc; clear all; close all

channelWidth = 4;       % channel width 1, 2, 4, 8, 16 MHz
param = 6;              % SISO: 0: BPSK, 1,2: QPSK, 3,4: 16QAM, 5,6,7: 64QAM, 8,9: 256QAM, 10: BPSK 2x repetition; 
                        % MIMO: 0: BPSK, 1,2: QPSK
GI = 'normal';          % guard interval SISO: normal, short; MIMO: normal
model = 'SISO';         % model type: SISO, MIMO

n = 10000;              % number of input bits

CNR = 30;               % Carrier-to-Noise ratio in units dB
% CNR = 0:1:34;     
channelType = 'AWGN';   % AWGN, Rician , Rayleigh
useEkvaliz = '0';       % 0 - Equalization OFF ; 1 - Equalization ON 

err = 1e-6;             % error floor to plot BER graphs 
end

fc = 863000000;         % carrier frequency in Hz
fs = 31250;             % sampling frequency in Hz

%% The call of the mcss.m file  

[Nsd, Nsa, Nbpsc, Ncbps, Ndbps, Kmod, pilots] = mcss(channelWidth, param, model); % modulation and coding schemes (MCSS)

%% To generate data bits - Data Input

PSDU = randi([0 1], 1, n);                       % generate binary random matrix

DATA = [zeros(1, 16) PSDU zeros(1, 6)];          % adding 16 service bits and 6 tail bits

dataLength = length(DATA);
totalDataNum = ceil(dataLength / Ndbps) * Ndbps; % number of bits in data field
nPadding = totalDataNum - dataLength;            % number of pad bits
DATA = [DATA zeros(1, nPadding)];                % adding pad bits to data field

clear dataLength totalDataNum 

%% Scrambling

% data scrambling - the calling of the scrambler.m file
scrambledData = scrambler(DATA);

% replacing 6 scrambled tail bits by 6 zero bits
scrambledData(end-nPadding-5:end-nPadding) = zeros(1, 6);

%% Convolutional encoder

K = 7;                  % Constraint length of the code
codeGen = [133 171];    % Code generator for the convolutional encoder

% convert convolutional code polynomials to trellis description
trellis = poly2trellis(K, codeGen); 

% encodes the binary vector "scrambledData" using the convolutional encoder whose MATLAB trellis structure is "trellis"
switch (param)
    case {0, 1, 3}
        codedData = convenc(scrambledData, trellis); 
    case {2, 4, 6, 8}
        puncpat = [1 1 1 0 0 1];            % puncture pattern 3/4
        codedData = convenc(scrambledData, trellis, puncpat); 
    case 5
        puncpat = [1 1 1 0];                % puncture pattern 2/3
        codedData = convenc(scrambledData, trellis, puncpat); 
    case {7, 9}
        puncpat = [1 1 1 0 0 1 1 0 0 1];    % puncture pattern 5/6
        codedData = convenc(scrambledData, trellis, puncpat); 
    case 10
        dataReady = convenc(scrambledData, trellis);
        Cout = reshape(dataReady, 12, []);  % data must be repeated after each 12 bits
        sequenceS = repmat([1 0 0 0 0 1 0 1 0 1 1 1]', 1, size(Cout, 2));
        dataXor = xor(Cout, sequenceS);     % repeated data must be xor
        dataDouble = [Cout; dataXor];       % coded and repeated data
        codedData = reshape(dataDouble, 1, []);
    otherwise
        error(['Unknown mode ' num2str(param)]);
end

clear Cout sequenceS dataXor dataDouble
%% Interleaving

nCodedBits = length(codedData);                 % number of coded bits after convolution coding
switch channelWidth
    case {1, 2, 4, 8}
        x = reshape(codedData, Ncbps, []);      % spolecne se proklada vzdy jen Ncbps bitu
    case 16
        x = reshape(codedData, Ncbps/2, []);    % for 16 MHz only Ncbps/2 bits
end
[nRow, nCol] = size(x);

firstIntrlvd = zeros(nRow, nCol);
secIntrlvd = zeros(nRow, nCol);

switch channelWidth
    case 1
        Ncol = 8;       % the width of the interleaver
    case 2
        Ncol = 13;
    case 4
        Ncol = 18;
    case {8, 16}
        Ncol = 26;
end
Nrows = nRow/Ncol;      % number of rows for OFDM symbol 

  
for i = 1:nCol
    firstIntrlvd(:,i) = matintrlv(x(:,i), Nrows, Ncol);                  % Reorder symbols by filling matrix by rows and emptying it by columns
    secIntrlvd(:,i) = longRunsIntrlv(firstIntrlvd(:,i), Nbpsc, Ncol);    % second permutation - the call of a custom function
end  

intrlvdData = reshape(secIntrlvd, 1, []);

clear i x firstIntrlvd secIntrlvd

%% Subcarrier Modulation Mapping (Constellation Mapping)

dataMatrix = reshape(intrlvdData, [Nbpsc,(nCodedBits)/Nbpsc])';   % reshape into block of Nbpsc colums 
decadicData = bi2de(dataMatrix, 'left-msb');                      % convert binary vectors according to Nbpsc to decimal numbers % convert binary vectors to decimal numbers
d = reshape(decadicData, Nsd, []);
modOrder = 2^Nbpsc;                                               % modulation order

% map on the IQ plane
switch (param)
    case {0, 10}
        dataMapped = pskmod(d, modOrder, 'gray') * Kmod * (-1);
    case {1, 2, 3, 4, 5, 6, 7, 8, 9}
        dataMapped = conj(qammod(d, modOrder, 'gray')) * Kmod;        
    otherwise
        error(['Unknown mode ' num2str(param)]);
end

if ~(exist('usedGUI','var')) && length(CNR) == 1
    plotIQdiag( dataMapped, modOrder )                             % plot constellation diagram (TX)
end

clear dataMatrix decadicData d
%% Space-Time Block Coding (STBC)

switch model
    case 'MIMO'
        Nss = 2;
        resMap = reshape(dataMapped, 2, []);       
        STBCMatrix(1,:,Nss-1) = resMap(1,:);
        STBCMatrix(1,:,Nss) = resMap(2,:);
        STBCMatrix(2,:,Nss-1) = -conj(resMap(2,:));
        STBCMatrix(2,:,Nss) = conj(resMap(1,:));        
        dataForOFDM = reshape(STBCMatrix, Nsd, [], 2);
    case 'SISO'
        Nss = 1;
        dataForOFDM = dataMapped;
end

clear resMap STBCMatrix

%% OFDM Modulation and Insertion of Pilots

pilots = pilots + (Nsa/2+1);
nPilot = length(pilots);

Pseq = [...       % Pilot sequence of 127 BPSK symbols
  1,1,1,1, -1,-1,-1,1, -1,-1,-1,-1, 1,1,-1,1, -1,-1,1,1, -1,1,1,-1,...
  1,1,1,1, 1,1,-1,1, 1,1,-1,1, 1,-1,-1,1, 1,1,-1,1, -1,-1,-1,1,...
  -1,1,-1,-1, 1,-1,-1,1, 1,1,1,1, -1,-1,1,1, -1,-1,1,-1, 1,-1,1,1,...
  -1,-1,-1,1, 1,-1,-1,-1, -1,1,-1,-1, 1,-1,1,1, 1,1,-1,1, -1,1,-1,1,...
  -1,-1,-1,-1, -1,1,-1,1, 1,-1,1,-1, 1,1,1,-1, -1,1,-1,-1, -1,1,1,1,...
  -1,-1,-1,-1, -1,-1,-1
  ];

[nx, nOFDM, ns] = size(dataForOFDM);
copies = ceil(nOFDM/length(Pseq));      % Number of copies needed
nPilotCop = repmat(Pseq, 1, copies);    % Generate enough copies

% the creation of the OFDM symbols
switch channelWidth
    case {1, 2}
        zerosDC = 1;
    case {4, 8}
        zerosDC = 3;
    case 16
        zerosDC = 17;
end

zero = Nsa - (Nsd+nPilot+zerosDC);      % to compute the number of zeros at the end of OFDM symbols
a = ceil(zero/2);                       % the number of zeros for the first half of an OFDM symbol
pPrototyp = [1 1 1 -1 -1 1 1 1 1 1 1 -1 -1 1 1 1];
pilotMatrix = zeros(Nsa, nOFDM, Nss);
allPilots = zeros(nPilot, nOFDM, Nss);

for s = 1:Nss
    addZerosOFDM(:,:,s) = [zeros(a,nOFDM,1); dataForOFDM(:,:,s); zeros(zero-a+nPilot+zerosDC,nOFDM,1)];
end
dataForOFDM = addZerosOFDM;
    
for ss = 1:Nss
    for i = 1:nPilot
        allPilots(i,:,ss) = nPilotCop(1:nOFDM) * pPrototyp(i);
        pilotMatrix(pilots(i),:,ss) = allPilots(i,:,ss);

        dataForOFDM(pilots(i)+1:end,:,ss) = dataForOFDM(pilots(i):end-1,:,ss);
        dataForOFDM(pilots(i),:,ss) = 0;

        switch channelWidth
            case {1, 2}
                if i == nPilot/2
                    dataForOFDM(Nsa/2+2:end,:,ss) = dataForOFDM(Nsa/2+1:end-1,:,ss);
                    dataForOFDM(Nsa/2+1,:,ss) = 0;
                end
            case {4, 8}
                if i == nPilot/2
                    dataForOFDM(Nsa/2+3:end,:,ss) = dataForOFDM(Nsa/2:end-3,:,ss);
                    dataForOFDM(Nsa/2:Nsa/2+2,:,ss) = 0;
                end            
            case 16
                if i == 4
                    dataForOFDM(130:end,:,ss) = dataForOFDM(127:end-3,:,ss);
                    dataForOFDM(127:129,:,ss) = 0;   
                elseif i == 8
                    dataForOFDM(262:end,:,ss) = dataForOFDM(251:end-11,:,ss);
                    dataForOFDM(251:261,:,ss) = 0;   
                elseif i == 12
                    dataForOFDM(387:end,:,ss) = dataForOFDM(384:end-3,:,ss);
                    dataForOFDM(384:386,:,ss) = 0;                   
                end
        end
    end
    
    ofdmFrame(:,:,ss) = dataForOFDM(:,:,ss) + pilotMatrix(:,:,ss);
    ofdmFrame(:,:,ss) = [ofdmFrame(Nsa/2+1:end,:,ss); ofdmFrame(1:Nsa/2,:,ss)];    % swap half OFDM symbol
    
end

%% IFFT and CP Insertion

% to provide Inverse Fast Fourier Transform (coversion from the frequency to the time domain)
ifftDataPar = ifft(ofdmFrame); 

% cyklic prefix (CP)
switch model
    case 'SISO'
        switch channelWidth
            case 1
                switch GI
                    case 'normal'
                        cycPrefix = ifftDataPar((end-floor(Nsd/4)+1):end,:,:);
                    case 'short'
                        cycPrefix = ifftDataPar((end-floor(Nsd/12)+1):end,:,:);
                    otherwise
                        error(['Unknown guard interval ' num2str(GI) ' for SISO']);
                end
            case {2, 4, 8, 16}
                switch GI
                    case 'normal'
                        cycPrefix = ifftDataPar((end-floor(Nsd/4)+1):end,:,:); 
                    otherwise   
                        error(['Unknown guard interval ' num2str(GI) ' for SISO']);
                end
        end   
case 'MIMO'
    switch GI
        case 'normal'
            cycPrefix = ifftDataPar((end-floor(Nsd/4)+1):end,:,:);
        otherwise
            error(['Unknown guard interval ' num2str(GI) ' for MIMO']);
    end
end

ifftDataPar = [cycPrefix; ifftDataPar];         % add CP

ifftDataSer = reshape(ifftDataPar, 1, [], Nss); % serial-to-parallel data

% Measurment of the signal power (before upsampling)
lFilt = 10*log10(0.8);    
signalPower = 10*log10((sum(ifftDataSer(:).*conj(ifftDataSer(:))))/length(ifftDataSer(:)));
signalPower = signalPower + 10*log((Nsd-nPilot)/Nsa) + 2*lFilt;   % ratio of used carriers vs all carriers

clear addZerosOFDM s ss pilotMatrix nPilotCop pilotMatrix copies Pseq pPrototyp lFilt ns

%% I/Q Modulator and Data Transmission

% real and imaginary part
dataI = real(ifftDataSer);
dataQ = imag(ifftDataSer);

% to increase sampling rate by upsampling
upsampling = 10;
dataII = upsample(dataI, upsampling);
dataQQ = upsample(dataQ, upsampling);

% raised cosine FIR pulse-shaping filter (RCF filter)
rollof  = 0.13;
span    = 60;
sps     = upsampling;
b = rcosdesign(rollof, span, sps, 'sqrt');

% to add zeros for filter delay
dataII(1,end+span*sps/2,:) = 0; 
dataQQ(1,end+span*sps/2,:) = 0;

% data filtering
datamodI = filter(b, 1, dataII); 
datamodQ = filter(b, 1, dataQQ);

% compensation of the filter delay
datamodII = datamodI(1,(sps/2)*span+1 : end,:);
datamodQQ = datamodQ(1,(sps/2)*span+1 : end,:);

% To cretate a carrier sginal
deltaf = fs * upsampling; % 31.25 kHz
t = 0 : (1/deltaf) : (length(datamodII)-1)/deltaf;
carrierI = cos(2*pi*t*fc);
carrierQ = sin(2*pi*t*fc);

% TX signal and transmission
signalTX = zeros(1,size(carrierI,2),Nss);
for cc = 1:Nss
    signalTX(:,:,cc) = datamodII(:,:,cc).*carrierI + datamodQQ(:,:,cc).*carrierQ;
end

% if ~(exist('usedGUI','var'))
%     figure
%     pwelch(signalTX(:,:,1));
% end

clear dataI dataQ dataII dataQQ datamodI datamodQ datamodII datamodQQ deltaf t

%% Channel
% the calling of channel models

if ~strcmp(channelType, 'AWGN')
    [signalTXChan] = channel( signalTX, model, channelType);
    signalTX = signalTXChan;

%     figure
%     pwelch(signalTXChan);
end

% MIMO
if strcmp(model, 'MIMO') && strcmp(channelType, 'AWGN')
    h0=1;h1=1;h2=1;h3=1;
	signalTX(:,:,1) = h0*signalTX(:,:,1) + h1*signalTX(:,:,2); 
	signalTX(:,:,2) = h2*signalTX(:,:,1) + h3*signalTX(:,:,2);
end


%% AWGN channel

BER = 1;

for nn = 1:length(CNR)
if BER > err 

SNRVec = CNR(nn);

if (exist('usedGUI','var'))
    set(handles.textStatus, 'String', ['Výpoèet BER pro SNR = ' num2str(SNRVec) 'dB']);
    drawnow
else
    disp(['Computing BER for SNR = ' num2str(SNRVec) 'dB'])
end

% add AWGN
signalRX = zeros(1,size(carrierI,2),Nss);
for aa = 1:Nss
    signalRX(:,:,aa) = awgn(signalTX(:,:,aa),SNRVec,signalPower);
end

%% Receiver (RX)

%% I/Q demodulator

for aa = 1:Nss
    dataRXI(:,:,aa) = signalRX(:,:,aa).*(carrierI);
    dataRXQ(:,:,aa) = signalRX(:,:,aa).*(carrierQ);
end

% raised cosine FIR pulse-shaping filter
b = rcosdesign(rollof, span, sps, 'sqrt');

% to add zeros for filter delay
dataRXI(1,end+span*(sps/2),:) = 0;
dataRXQ(1,end+span*(sps/2),:) = 0;

dataRXI = filter(b, 0.5, dataRXI);
dataRXQ = filter(b, 0.5, dataRXQ);

% compensation of filter delay
dataRXII = dataRXI(1,(sps/2)*span+1 : end,:);
dataRXQQ = dataRXQ(1,(sps/2)*span+1 : end,:);

% downsampling
dataRXIII = downsample(dataRXII, upsampling);
dataRXQQQ = downsample(dataRXQQ, upsampling);

for cc = 1:Nss
dataRX(:,:,cc) = dataRXIII(:,:,cc) + 1j*dataRXQQQ(:,:,cc);
end

clear aa dataRXI dataRXQ dataRXII dataRXQQ dataRXIII dataRXQQQ

%% Remove CP and to provide FFT

ifftDataParRecCyc = reshape(dataRX, Nsa+length(cycPrefix(:,1)), [], Nss);   % parallel-to-serial data
ifftDataParRec = ifftDataParRecCyc(length(cycPrefix(:,1))+1:end,:,:);       % remove CP

fftData = fft(ifftDataParRec);                                              % to provide FFT
fftData = [fftData(Nsa/2+1:end,:,:); fftData(1:Nsa/2,:,:)];

% to plot OFDM spectrum - to hide it, you must commented the code from "if" to "end"
if ~(exist('usedGUI','var')) && length(CNR) == 1
    BBDataCH = reshape(ifft(fftData),1,[],1);
    N = 8192;               % choose window for plot
    while N > length(BBDataCH)
        N = N/2;
    end
    
    figure
    Pxx = abs(fft(BBDataCH(1,1:N)))/N;
    df = fs/(N/Nsa)/1e6;
    xf = -fs*Nsa/2e6:df:fs*Nsa/2e6-df;
    plot(xf,10*log(Pxx/max(Pxx)));
    xlabel('Frequency [MHz]')
    ylabel('Normalized magnitude [dB]')
end

fftDataa = fftData;

%% Equalization (Zero-Force)
% the call of equalizer.m file

if strcmp(useEkvaliz, '1') && strcmp(model, 'SISO')
	ekvalizChan = equalizer(fftData, pilots, allPilots, Nsa);
	fftDataa = ekvalizChan;
end

    i = nPilot;
    while i >= 1    
        fftDataa(pilots(i):end-1,:,:) = fftDataa(pilots(i)+1:end,:,:);

        switch channelWidth
            case {1, 2}
                if i == (nPilot/2)+1
                    fftDataa(Nsa/2+1:end-1,:,:) = fftDataa(Nsa/2+2:end,:,:);
                end
            case {4, 8}
                if i == (nPilot/2)+1
                    fftDataa(Nsa/2:end-3,:,:) = fftDataa(Nsa/2+3:end,:,:);
                end
            case 16
                if i == 5
                    fftDataa(127:end-3,:,:) = fftDataa(130:end,:,:);  
                elseif i == 9
                    fftDataa(251:end-11,:,:) = fftDataa(262:end,:,:);
                elseif i == 13
                    fftDataa(384:end-3,:,:) = fftDataa(387:end,:,:);                 
                end
        end

        i = i-1;
    end
    
fftDataaa = fftDataa(a+1:end-(zero-a+nPilot+zerosDC),:,:);

spatialMapRec = fftDataaa;

clear ifftDataParRecCyc ifftDataParRec fftDataa i N BBDataCH

%% STBC decoding

switch model
    case 'MIMO'
        if ~(exist('usedGUI','var')) && length(CNR) == 1     
            plotIQdiag(spatialMapRec(:,:,1), modOrder) % plot constellation recieve data
        end
        % recieved data for antennas Rx1, Rx2 as rows
        resAL = reshape(spatialMapRec,2,[],2);
        h0=1;h1=1;h2=1;h3=1;

        for n=1:length(resAL);
            AL(1,n) = (conj(h0)*resAL(1,n,1))+(h1*conj(resAL(2,n,1)))+(conj(h2)*resAL(1,n,2))+(h3*conj(resAL(2,n,2)));           
            AL(2,n) = (conj(h1)*resAL(1,n,1))-(h0*conj(resAL(2,n,1)))+(conj(h3)*resAL(1,n,2))-(h2*conj(resAL(2,n,2)));
        end;
        
        norm = (abs(h0)^2 + abs(h1)^2 + abs(h2)^2 + abs(h3)^2);
        AL = AL./norm;
        dataMappedRec = reshape(AL,Nsd,[]);

    case 'SISO'
        dataMappedRec = spatialMapRec;
end

% to compute MER
if strcmp(model, 'SISO')
    errorVector = (dataMappedRec - dataMapped);
    errorRMS = sum(abs(errorVector).^2)/length(errorVector);
    RMS = sum(abs(dataMapped).^2)/length(dataMapped);
    MERdB(nn) = 10*log10(RMS/errorRMS);
end

clear AL resAL errorVector errorRMS RMS

%% Constellation Demapping

if ~(exist('usedGUI','var')) && length(CNR) == 1
    plotIQdiag(dataMappedRec, modOrder)                         % plot constellation recieve data
end

switch (param)
    case {0, 10}
        z = pskdemod(dataMappedRec./Kmod*(-1), 2^Nbpsc, 'gray');       
    case {1, 2, 3, 4, 5, 6, 7, 8, 9}
        z = qamdemod(conj(dataMappedRec./Kmod), 2^Nbpsc, 'gray');       
    otherwise
        error(['Unknown mode ' num2str(param)]);
end

b = de2bi(z, 'left-msb');                                       % convert decimal numbers to binary row vectors
dataVector = reshape(b', [1,nCodedBits]);                       % reshape into one row vector

%% Deinterleaving

switch channelWidth
    case {1, 2, 4, 8}
        xx = reshape(dataVector, Ncbps, []);
    case 16
        xx = reshape(dataVector, Ncbps/2, []);
end

        secDeintrlvd = zeros(nRow, nCol);
        firstDeintrlvd = zeros(nRow, nCol);
        for i = 1:nCol
            secDeintrlvd(:,i) = longRunsDeintrlv(xx(:,i), Nbpsc, Ncol);    % the call of a custom function
            firstDeintrlvd(:,i) = matdeintrlv(secDeintrlvd(:,i), Nrows, Ncol);
        end

deintrlvdData = reshape(firstDeintrlvd, 1, []);

clear xx 

%% Viterbi Decoder

tblen = 30;     % traceback length 

%decodes the vector using the Viterbi algorithm
switch param
    case {0, 1, 3}
        [numErrs, BER(nn)] = biterr(deintrlvdData, codedData);
        decoded = vitdec(deintrlvdData,trellis,tblen,'trunc','hard');
    case {2, 4, 5, 6, 7, 8, 9}
        [numErrs, BER(nn)] = biterr(deintrlvdData, codedData);
        decoded = vitdec(deintrlvdData,trellis,tblen,'trunc','hard', puncpat);
    case 10
        dataDouble = reshape(deintrlvdData, 24, []);
        Cout = dataDouble(1:12,:);
        CoutLine = reshape(Cout, 1, []);
        [numErrs, BER(nn)] = biterr(CoutLine, dataReady);
        decoded = vitdec(CoutLine,trellis,tblen,'trunc','hard');
end

clear dataDouble Cout CoutLine

%% Descrambling
% the call of function descrambler.m 

plainData = descrambler(decoded);

%% Calcualtion of BER before and after Viterbi decoding and show MER

[numErrss, BERafterViterbi(nn)] = biterr(plainData(1:(end-nPadding-6)),DATA(1:(end-nPadding-6)));

end
end

% to show the values of BER and MER
if (exist('usedGUI','var'))
    set(handles.textStatus, 'String', 'Simulace dokonèena');
    if length(BER) == 1
        set(handles.BERbeforeViterbi,'String', num2str(BER));
        set(handles.BERafterViterbi,'String', num2str(BERafterViterbi));
        switch model
            case 'SISO'
                set(handles.textMER,'String', num2str(MERdB));
            case 'MIMO'
                set(handles.textMER,'String', '');
        end
    else
        set(handles.BERbeforeViterbi,'String', num2str(BER(end)));
        set(handles.BERafterViterbi,'String', num2str(BERafterViterbi(end)));
        switch model
            case 'SISO'
                set(handles.textMER,'String', num2str(MERdB(end)));
            case 'MIMO'
                set(handles.textMER,'String', '');
        end
    end
    drawnow
else
    disp('Simulation finished')
    if length(CNR) > 1
        figure
        grid on
        semilogy(CNR(1:length(BER)),BER)  
        hold on
        semilogy(CNR(1:length(BER)),BERafterViterbi)  
    elseif length(BER) == 1
        BER
        BERafterViterbi
        if strcmp(model, 'SISO')
            MERdB
        end
    end
end
