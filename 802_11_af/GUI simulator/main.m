function varargout = main(varargin)
% MAIN MATLAB code for main.fig
%      MAIN, by itself, creates a new MAIN or raises the existing
%      singleton*.
%
%      H = MAIN returns the handle to a new MAIN or the handle to
%      the existing singleton*.
%
%      MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN.M with the given input arguments.
%
%      MAIN('Property','Value',...) creates a new MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main

% Last Modified by GUIDE v2.5 10-Jun-2020 15:32:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_OpeningFcn, ...
                   'gui_OutputFcn',  @main_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before main is made visible.
function main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main (see VARARGIN)
% Choose default command line output for main
handles.output = hObject;

%axes(handles.axes3)
%matlabImage = imread('transmitter.png');
%matlabImage = imread('80211af.png');
%image(matlabImage)
%axis off
%axis image

%slider setup
%set(handles.slider1, 'min', 54); % slider minimum
%set(handles.slider1, 'max', 1000); % slider maximum
%set(handles.slider1, 'Value', 54); % slider actual value
%minSliderValue = get(handles.slider1, 'Min');
%maxSliderValue = get(handles.slider1, 'Max');
%theRange = maxSliderValue - minSliderValue;
%steps = [1/theRange, 10/theRange];
%set(handles.slider1, 'SliderStep', steps);

value_radiobutton = get(handles.radiobutton8,'Value');
if (value_radiobutton == 1)
    set(handles.checkbox1, 'Enable', 'off');
else
    set(handles.checkbox1, 'Enable', 'on');
end

% Update handles structure
guidata(hObject, handles);

updateParameters(hObject,eventdata, handles)

% UIWAIT makes main wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function input1_Callback(hObject, eventdata, handles)
% hObject    handle to input1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input1 as text
%        str2double(get(hObject,'String')) returns contents of input1 as a double


% --- Executes during object creation, after setting all properties.
function input1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% prepation for simulation
clc;
if(ishandle(1))
    hFig = figure(1);
    close(hFig)
end
hMain = findobj(0, 'type', 'figure');
movegui(hMain,'northwest');

% set up timer and status bar
tic % start timing
perc = 0;
w = waitbar(perc/100,sprintf('%d%% Initialization parameters...',perc)); % initialize waitbar
% set(w,'WindowStyle', 'modal');
% frames = java.awt.Frame.getFrames();
% frames(end).setAlwaysOnTop(1);

%% Read GUI user inputs

% gets selected antenna system cofiguration by user from popupmenu3
configurations = get(handles.popupmenu3,'string');
selected_configuration = get(handles.popupmenu3,'Value');
antenna_mode = configurations{selected_configuration};

switch antenna_mode
    case '1x1 SISO' % User selects 1x1 SISO
        N_Tx = 1; % Number of Tx antennas
        N_Rx = 1; % Number of spatial streams
    case '2x1 MISO' % User selects 2x1 MISO
        N_Tx = 2; % Number of Tx antennas
        N_Rx = 1; % Number of Rx antennas
    case '2x2 MIMO' % User selected 2x2 MIMO
        N_Tx = 2; % Number of Tx antennas 
        N_Rx = 2; % Number of Rx antennas
end

N_SS = max([N_Tx N_Rx]); % Number of spatial streams

% gets selected BCU bandwidth by user from popupmenu4
BCU_bandwidths = get(handles.popupmenu4,'string');
selected_BCU_bandwidth = get(handles.popupmenu4,'Value');
BCU_bandwidth = BCU_bandwidths{selected_BCU_bandwidth};

% gets selected channel configuration by user from popupmenu5
CH_bandwidths = get(handles.popupmenu5,'string');
selected_CH_bandwidth = get(handles.popupmenu5,'Value');
CH_bandwidth = CH_bandwidths{selected_CH_bandwidth};

% gets transmission format choice defined in uibuttongroup7
selected_radiobutton = get(handles.uibuttongroup7,'SelectedObject');
selected_format = get(selected_radiobutton,'string');
result = strcmp(selected_format,'Non-HT');
if (result == 1)
    transmission_type = 'non-HT';
    
    % read Modulation Coding Scheme parameters for duplicate non-HT mode
    T=readtable('MCS.txt','Format','%s %s %s %s %f %f %f %s');
else
    transmission_type = 'VHT';
    
    % read Modulation Coding Scheme parameters for duplicate VHT mode
    T=readtable('MCS_TVHT.txt','Format','%s %s %s %f %f %f %s');
end

T.Properties.RowNames = T.MCS'; % define MCS rows

% gets selected MCS index by user from popupmenu1
MCS_indexes = get(handles.popupmenu1,'string');
selected_MCS_index = get(handles.popupmenu1,'Value');
MCS_index = MCS_indexes{selected_MCS_index};

% gets number of data bytes defined by user from input1
number_of_data_bytes = get(handles.input1,'string');
[num,status]=str2num(number_of_data_bytes); % check if input string is number
if (status == 1) % data_input is correct
    
% gets cyclic prefix choice defined in uibuttongroup4
selected_radiobutton = get(handles.uibuttongroup4,'SelectedObject');
selected_GI = get(selected_radiobutton,'string');
result = strcmp(selected_GI,'Normal');
if (result == 1)
    GI = 'GI';
else
    GI = 'GIS';
end

% gets SNR value in uibuttongroup5
selected_radiobutton = get(handles.uibuttongroup5,'SelectedObject');
selected_SNR = get(selected_radiobutton,'string');
result = strcmp(selected_SNR,'Static');

if (result == 1) % if SNR static
    
    SNR_value = get(handles.edit7,'string'); % gets SNR value from edit7
    [SNR,status1]=str2num(SNR_value); % check if input string is number
    if (status1 == 0)
        error('SNR value must be number!');
    end
    
else % if SNR dynamic
    
    SNR_start_input = get(handles.edit4,'string'); % gets SNR value from edit4
    SNR_stop_input = get(handles.edit5,'string'); % gets SNR value from edit5
    SNR_step_input = get(handles.edit6,'string'); % gets SNR value from edit6
    
    [SNR_start,status2]=str2num(SNR_start_input); % check if input string is number
    [SNR_stop,status3]=str2num(SNR_stop_input); % check if input string is number
    [SNR_step,status4]=str2num(SNR_step_input); % check if input string is number
    if ( status2 == 0 || status3 == 0 || status4 == 0 )
        error('SNR value must be number!');
    end   
    
    % SNR vector
    SNR = SNR_start:SNR_step:SNR_stop;
end
% gets RF channel model in uibuttongroup6
selected_radiobutton = get(handles.uibuttongroup6,'SelectedObject');
selected_RF_model = get(selected_radiobutton,'string');

% gets equalization choice in checkbox1
if (get(handles.checkbox1, 'Value'))
    equalizer_flag = 1;
else
    equalizer_flag = 0;
end

%% Variables initialization

    BCU_W = handles.BCU_W; % one BCU size
    FFT_size = handles.IFFT_size; % size of IFFT/FFT window
    f_BW = handles.BW * (10 ^ 6); % Bandwidth 8 MHz
    delta_f = f_BW/BCU_W; % Subcarrier frequency spacing 0.0625 MHz
    T_FFT = 1/delta_f; % IFFT/FFT period [micros]
    N_SD = handles.N_SD; % Number of data subcarriers per BCU 6/7/8 MHz
    N_SP = handles.N_SP; % Number of pilot subcarriers
    N_ST = handles.N_ST; % Number of subcarriers total

    switch(GI)
        case 'GI'
            N_GI = floor(BCU_W/4); % GI duration
            T_GI = T_FFT/4; % OFDM symbol period 0.8 [micros]
        case 'GIS' 
            N_GI = floor(BCU_W/8); % GI short duration
            T_GI = T_FFT/8; % OFDM symbol period 0.4 [micros]
        otherwise error('Incorrect guard interval! Could be GI/GIS')
    end

    T_OFDM = T_FFT + T_GI; % full OFDM symbol period [micros]
    fs = 1/T_OFDM; % OFDM symbol frequency [Hz]
    fs = delta_f; % symbol frequency spacing [Hz]
    N_OFDM = BCU_W + N_GI; % OFDM symbol length
    s = 127; % scrambler/descrambler initialization "1111111"

    
    fc = 868e+6; % Carrier frequency [Hz]
    
    % filter parameters
    %delay = 100; % RRC filter time delay
    beta = 0.03; % RRC filter roll off factor
    span = 60; % filter span
    upsampling = 10; % oversampling factor
    downsampling = upsampling; % downsampling coefficient

perc = 10;
waitbar(perc/100,w,sprintf('%d%% Transmitter processing...',perc));

%% Generate SIGNAL + DATA parts of PPDU (IEEE 802.11af)

    PSDU=PSDU_generator(num); % generates random user data bytes PSDU with legth of number_of_data_bytes      
    
    if ( strcmp(transmission_type,'non-HT') )
        SIGNAL=PPDU_SIGNAL(MCS_index,num); % generate SIGNAL field using data_rate and number of user bytes
        DATA=PPDU_DATA(PSDU,SIGNAL,T); % generate DATA field
    elseif ( strcmp(transmission_type,'VHT') )
        DATA=PPDU_DATA_TVHT(PSDU, MCS_index, T); % generate DATA field       
    end   
    
%% Scrambler DATA

    [scrambled_data, half_byte] = scrambler(DATA,s); % scrambles DATA input
    
%% FEC encoding + puncturing DATA + SIGNAL

    % check MCS vector in table
    T.Properties.RowNames = T.MCS'; % define MCS index rows
    R = T(MCS_index,:).R; %read coding rate from T
    R = R{1}; % converts cell to string

%     if ( strcmp(transmission_type,'non-HT') )
%         R = T(selected_MCS,:).R; %read coding rate from T
%         R = R{1}; % converts cell to string
%     elseif ( strcmp(transmission_type,'VHT') )
%         MCS = MCS + 1; % for correct reading MCS index from table
%         R = T(selected_MCS,:).R; %read coding rate from T
%         R = R{1}; % converts cell to string
%     end
    
    coded_data = fec_coder(scrambled_data,R); % convolutionally coded data input 
    
    if ( strcmp(transmission_type,'non-HT') )
        coded_signal = fec_coder(SIGNAL,'1/2'); % convolutionally coded signal input  
    end
    
%    tmp = stream_parser_TVHT(coded_data);
    
if ( ( N_Tx == 1 ) && ( N_Rx == 1 ) ) % if used SISO mode
    
%% Interleaving DATA + SIGNAL

    % check MCS vector in table
    T.Properties.RowNames = T.MCS'; % define MCS index rows
    Modulation = T(MCS_index,:).Modulation; %read coding rate from T
    Modulation = Modulation{1}; % converts cell to string

%     if ( strcmp(transmission_type,'non-HT') )
%         Modulation = T(data_rate,:).Modulation; % read modulation from T
%         Modulation = Modulation{1}; % converts cell to string
%     elseif ( strcmp(transmission_type,'VHT') )
%         Modulation = T2(MCS,:).Modulation; % read modulation from T2
%         Modulation = Modulation{1}; % converts cell to string
%     end
    
    [idata1, idata2] = interleaver(coded_data, CH_bandwidth, transmission_type, Modulation); % interleaved data input 
    % idata1 - after 1st permutation | idata2 - after 2nd permutation
    
    if ( strcmp(transmission_type,'non-HT') )
        [idata3, idata4] = interleaver(coded_signal, CH_bandwidth, transmission_type, 'BPSK'); % interleaved signal input 
    end
%% Subcarrier modulation DATA + SIGNAL

    modulated_data = modulator(idata2,Modulation); %modulated data input (works for BPSK modulation)
    
    if ( strcmp(transmission_type,'non-HT') )
        modulated_signal = modulator(idata4,'BPSK'); %modulated signal input (works for BPSK modulation)
        modulated_output = [modulated_signal modulated_data]; % modulated DATA + SIGNAL output
    elseif ( strcmp(transmission_type,'VHT') )
        modulated_output = modulated_data; % modulated DATA output
    end
   

%% OFDM modulation (pitot insertion + adding GI) DATA + SIGNAL

    [ofdm_data, data_flag] = ofdm_modulator2(modulated_data,transmission_type,CH_bandwidth,BCU_bandwidth,GI); % DATA OFDM modulation
    
    if ( strcmp(transmission_type,'non-HT') ) 
        [ofdm_signal, signal_flag] = ofdm_modulator2(modulated_signal,transmission_type,CH_bandwidth,BCU_bandwidth,GI); % SIGNAL OFDM modulation
         ofdm_output = [ofdm_signal ofdm_data]; % complete DATA + SIGNAL OFDM symbols
    elseif ( strcmp(transmission_type,'VHT') )
        ofdm_output = ofdm_data; % complete DATA OFDM symbols
    end
    

    
%% Signal power measurement before upsampling

    lFilt = 10*log10(0.9); % 0.8  VHT
    
    signal_power = 10*log10((sum(ofdm_output.*conj(ofdm_output)))/length(ofdm_output));
    %signal_power = signal_power + 10*log(N_ST/BCU_W) + 2*lFilt;   % ratio of used carriers vs all carriers
    signal_power = signal_power + 10*log((N_SD-N_SP)/BCU_W) + 2*lFilt;   % ratio of used carriers vs all carriers

%% Pulse shaping (RRC filter) + IQ modulation


    [I_filt, Q_filt] = RRC_filter(ofdm_output,upsampling,beta,span); % RRC filtering
    
    [output_signal, carrier_real, carrier_imag ] = IQ_modulator(I_filt, Q_filt,upsampling,fc,delta_f); % IQ modulation

perc = 40;
waitbar(perc/100,w,sprintf('%d%% Transmitter processing...',perc));
perc_step = floor(50/length(SNR));

%% RF fading channel models
   
 if (result == 1) % SNR static
     
      RF_signal = RF_channel(selected_RF_model,output_signal,SNR,signal_power,fs); 
 
%% IQ demodulator + Filter (RRC filter)

    [I,Q] = IQ_demodulator(RF_signal, carrier_real, carrier_imag); % IQ demodulation

    received_signal = RRC_filter2(I,Q,downsampling,beta,span); % RRC filtration

%% OFDM demodulator DATA + SIGNAL

if ( strcmp(transmission_type,'non-HT') ) 
    % split received signal into SIGNAL and DATA
    ofdm_signal_part = received_signal(1:N_OFDM); % ofdm SIGNAL
    ofdm_data_part = received_signal(N_OFDM+1:end); % ofdm DATA

    % BPSK modulated SIGNAL
    [ofdm_demodulated_signal, ofdm_demodulated_signal_pilots] = ofdm_demodulator2(ofdm_signal_part,transmission_type,CH_bandwidth,BCU_bandwidth,GI,signal_flag,equalizer_flag);
    % PSK/QAM modulated DATA 
    %equalizer_flag = 1;
    [ofdm_demodulated_data, ofdm_demodulated_data_pilots] = ofdm_demodulator2(ofdm_data_part,transmission_type,CH_bandwidth,BCU_bandwidth,GI,data_flag,equalizer_flag);
    %ofdm_demodulated_data = ofdm_demodulator2(ofdm_data_part,transmission_type,CH_bandwidth,BCU_bandwidth,GI,data_flag);
    
    ofdm_demodulated_output_pilots = [ofdm_demodulated_signal_pilots ofdm_demodulated_data_pilots]; % modulated SIGNAL + DATA
    ofdm_demodulated_output = [ofdm_demodulated_signal ofdm_demodulated_data]; % modulated SIGNAL + DATA
elseif ( strcmp(transmission_type,'VHT') ) 
    [ofdm_demodulated_data, ofdm_demodulated_data_pilots] = ofdm_demodulator2(received_signal,transmission_type,CH_bandwidth,BCU_bandwidth,GI,data_flag,equalizer_flag);
    ofdm_demodulated_output_pilots = ofdm_demodulated_data_pilots;
    ofdm_demodulated_output = ofdm_demodulated_data; % modulated DATA
end
%% Subcarrier demodulator DATA + SIGNAL

if ( strcmp(transmission_type,'non-HT') )
    demodulated_signal = demodulator(ofdm_demodulated_signal,'BPSK'); % interleaved SIGNAL
    demodulated_data = demodulator(ofdm_demodulated_data,Modulation); % interleaved DATA
elseif ( strcmp(transmission_type,'VHT') )
    demodulated_data = demodulator(ofdm_demodulated_data,Modulation); % interleaved DATA
end

%% Deinterleaving DATA + SIGNAL

if ( strcmp(transmission_type,'non-HT') )
    deinterleaved_signal = deinterleaver(demodulated_signal, CH_bandwidth, transmission_type, 'BPSK'); % deinterleaved SIGNAL
    deinterleaved_data = deinterleaver(demodulated_data, CH_bandwidth, transmission_type, Modulation); % deinterleaved DATA
elseif ( strcmp(transmission_type,'VHT') )
    deinterleaved_data = deinterleaver(demodulated_data, CH_bandwidth, transmission_type,Modulation); % deinterleaved DATA
end
%% FEC decoding + puncturing DATA + SIGNAL

if ( strcmp(transmission_type,'non-HT') )
    decoded_signal = fec_decoder(deinterleaved_signal,'1/2'); % convolutionally decoded signal
end
    decoded_data = fec_decoder(deinterleaved_data,R); % convolutionally decoded data

%% Descrambler DATA
    
    descramblered_data = descrambler(decoded_data,s,half_byte); % descrambler DATA

%% BER calculation

if ( strcmp(transmission_type,'non-HT') )    
    TxData = [coded_signal coded_data]; % transmitted bits
elseif ( strcmp(transmission_type,'VHT') )
    TxData = coded_data; % transmitted bits
end
    
    % pre Viterbi decoder
if ( strcmp(transmission_type,'non-HT') )
    RxData = [deinterleaved_signal deinterleaved_data]; % received bits before Viterbi decoder
elseif ( strcmp(transmission_type,'VHT') )
    RxData = deinterleaved_data; % received bits before Viterbi decoder
end
    nErrors = biterr(TxData,RxData); % calculate the number of bit errors
    numErrs = nErrors; % number of error bits
    numBits = length(TxData); % number of transmitted bits
    pre_BER_vector = numErrs / numBits % BER calculated

if ( strcmp(transmission_type,'non-HT') )
    TxData = [SIGNAL DATA]; % transmitted bits
elseif ( strcmp(transmission_type,'VHT') )
    TxData = DATA; % transmitted bits
end

    %tblen = 30;
    % post Viterbi decoder
if ( strcmp(transmission_type,'non-HT') )    
    RxData = [decoded_signal descramblered_data]; % received bits after Viterbi decoder
elseif ( strcmp(transmission_type,'VHT') )
    RxData = descramblered_data; % received bits after Viterbi decoder
end
    nErrors = biterr(TxData,RxData); % calculate the number of bit errors
    %nErrors = biterr(RxData(tblen+1:end),TxData(1:end-tblen));
    numErrs = nErrors; % number of error bits
    numBits = length(TxData); % number of transmitted bits
    post_BER_vector = numErrs / numBits % BER calculated
    
    perc = perc + perc_step;
    waitbar(perc/100,w,sprintf('%d%% Receiver processing...',perc));
    
%% MER + EVM calculation

    MER_vector = MER_calc(modulated_output, ofdm_demodulated_output); % MER calculation
    
    EVM_vector = EVM_calc(modulated_output, ofdm_demodulated_output); % EVM calculation
    
else % SNR dynamic
    
    pre_BER_vector = []; % vector for BER before Viterbi decoder on each SNR
    post_BER_vector = []; % vector for BER after Viterbi decoder on each SNR
    MER_vector = []; % vector for MER calculations on each SNR
    EVM_vector = []; % vector for EVM calculations on each SNR

    if ( strcmp(transmission_type,'non-HT') )    
        pre_TxData = [coded_signal coded_data]; % transmitted bits
    elseif ( strcmp(transmission_type,'VHT') )
        pre_TxData = coded_data; % transmitted bits
    end
    
    if ( strcmp(transmission_type,'non-HT') )
        post_TxData = [SIGNAL scrambled_data]; % transmitted bits
    elseif ( strcmp(transmission_type,'VHT') )
        post_TxData = scrambled_data; % transmitted bits
    end

    %pre_TxData = [coded_signal coded_data]; % transmitted bits
    %post_TxData = [SIGNAL scrambled_data]; % transmitted bits
    pre_numBits = length(pre_TxData); % number of transmitted bits ( 1x transmission )
    post_numBits = length(post_TxData); % number of transmitted bits ( 1x transmission )   
    
    for n=1:length(SNR)
        
        pre_RxData = 0;
        post_RxData = 0;
        
        RF_signal = RF_channel(selected_RF_model,output_signal,SNR(n),signal_power,fs);
        
        %% IQ demodulator + Filter (RRC filter)

        downsampling = upsampling; % downsampling coefficient

        [I,Q] = IQ_demodulator(RF_signal, carrier_real, carrier_imag); % IQ demodulation

        received_signal = RRC_filter2(I,Q,downsampling,beta,span); % RRC filtration

        %% OFDM demodulator DATA + SIGNAL
        
        if ( strcmp(transmission_type,'non-HT') ) 
            % split received signal into SIGNAL and DATA
            ofdm_signal_part = received_signal(1:N_OFDM); % ofdm SIGNAL
            ofdm_data_part = received_signal(N_OFDM+1:end); % ofdm DATA

            % BPSK modulated SIGNAL
            [ofdm_demodulated_signal, ofdm_demodulated_signal_pilots] = ofdm_demodulator2(ofdm_signal_part,transmission_type,CH_bandwidth,BCU_bandwidth,GI,signal_flag,equalizer_flag);
            % PSK/QAM modulated DATA 
            %equalizer_flag = 1;
            [ofdm_demodulated_data, ofdm_demodulated_data_pilots] = ofdm_demodulator2(ofdm_data_part,transmission_type,CH_bandwidth,BCU_bandwidth,GI,data_flag,equalizer_flag);
            %ofdm_demodulated_data = ofdm_demodulator2(ofdm_data_part,transmission_type,CH_bandwidth,BCU_bandwidth,GI,data_flag);

            ofdm_demodulated_output_pilots = [ofdm_demodulated_signal_pilots ofdm_demodulated_data_pilots]; % modulated SIGNAL + DATA
            ofdm_demodulated_output = [ofdm_demodulated_signal ofdm_demodulated_data]; % modulated SIGNAL + DATA
            
        elseif ( strcmp(transmission_type,'VHT') ) 
            [ofdm_demodulated_data, ofdm_demodulated_data_pilots] = ofdm_demodulator2(received_signal,transmission_type,CH_bandwidth,BCU_bandwidth,GI,data_flag,equalizer_flag);
            ofdm_demodulated_output_pilots = ofdm_demodulated_data_pilots;
            ofdm_demodulated_output = ofdm_demodulated_data; % modulated DATA
        end

        %% Subcarrier demodulator DATA + SIGNAL
        
        if ( strcmp(transmission_type,'non-HT') )
            demodulated_signal = demodulator(ofdm_demodulated_signal,'BPSK'); % interleaved SIGNAL
            demodulated_data = demodulator(ofdm_demodulated_data,Modulation); % interleaved DATA
        elseif ( strcmp(transmission_type,'VHT') )
            demodulated_data = demodulator(ofdm_demodulated_data,Modulation); % interleaved DATA
        end

        %% Deinterleaving DATA + SIGNAL

        %deinterleaved_signal = deinterleaver(demodulated_signal,'BPSK'); % deinterleaved SIGNAL
        %deinterleaved_data = deinterleaver(demodulated_data,Modulation); % deinterleaved DATA 
        
        if ( strcmp(transmission_type,'non-HT') )
            deinterleaved_signal = deinterleaver(demodulated_signal, CH_bandwidth, transmission_type, 'BPSK'); % deinterleaved SIGNAL
            deinterleaved_data = deinterleaver(demodulated_data, CH_bandwidth, transmission_type, Modulation); % deinterleaved DATA
        elseif ( strcmp(transmission_type,'VHT') )
            deinterleaved_data = deinterleaver(demodulated_data, CH_bandwidth, transmission_type, Modulation); % deinterleaved DATA
        end

        %% FEC decoding + puncturing DATA + SIGNAL

        if ( strcmp(transmission_type,'non-HT') )
            decoded_signal = fec_decoder(deinterleaved_signal,'1/2'); % convolutionally decoded signal
            decoded_data = fec_decoder(deinterleaved_data,R); % convolutionally decoded data
        elseif ( strcmp(transmission_type,'VHT') )
            decoded_data = fec_decoder(deinterleaved_data,R); % convolutionally decoded data
        end

        %% Descrambler DATA

        descramblered_data = descrambler(decoded_data,s,half_byte); % descrambler DATA
        
        %% BER calculation
        
        % pre Viterbi decoder
        %pre_RxData = [deinterleaved_signal deinterleaved_data]; % received bits before Viterbi decoder
        if ( strcmp(transmission_type,'non-HT') )
            pre_RxData = [deinterleaved_signal deinterleaved_data]; % received bits before Viterbi decoder
        elseif ( strcmp(transmission_type,'VHT') )
            pre_RxData = deinterleaved_data; % received bits before Viterbi decoder
        end
        nErrors = biterr(pre_TxData,pre_RxData); % calculate the number of bit errors
        numErrs = nErrors; % number of error bits
        %numBits = length(pre_TxData); % number of transmitted bits
        pre_BER_actual = numErrs / pre_numBits; % BER before Viterbi decoder calculated
        pre_BER_vector = [pre_BER_vector pre_BER_actual]; % fill BER vector with calculated BER before Viterbi decoder on each SNR
    
        % post Viterbi decoder
        %post_RxData = [decoded_signal decoded_data]; % received bits after Viterbi decoder
        if ( strcmp(transmission_type,'non-HT') )    
            post_RxData = [decoded_signal decoded_data]; % received bits after Viterbi decoder
        elseif ( strcmp(transmission_type,'VHT') )
            post_RxData = decoded_data; % received bits after Viterbi decoder
        end
        nErrors = biterr(post_TxData,post_RxData); % calculate the number of bit errors
        numErrs = nErrors; % number of error bits
        %numBits = length(TxData); % number of transmitted bits
        post_BER_actual = numErrs / post_numBits; % BER after Viterbi decoder calculated
        post_BER_vector = [post_BER_vector post_BER_actual]; % fill BER vector with calculated BER after Viterbi decoder on each SNR
        
        perc = perc + perc_step;
        waitbar(perc/100,w,sprintf('%d%% Receiver processing...',perc));
        
        %% MER + EVM calculation

        MER_act = MER_calc(modulated_output, ofdm_demodulated_output); % MER calculation
        MER_vector = [MER_vector MER_act]; % MER vector accumulation

        EVM_act = EVM_calc(modulated_output, ofdm_demodulated_output); % EVM calculation
        EVM_vector = [EVM_vector EVM_act]; % EVM vector accumulation
        
        
    end
    
 end

elseif ( ( N_Tx == 2 ) && ( N_Rx == 1 ) ) % if used MISO mode

    %% Stream parser 
    
    N_CBPS = T(data_rate,:).N_CBPS; % read number of coded bits per symbol from T
    coded_data_blocks = length(coded_data) / N_CBPS; % number of coded data blocks
    
    coded_Tx1 = coded_signal; % 1st coded Tx stream (SIGNAL + DATA) 
    coded_Tx2 = []; % 2nd coded Tx stream (DATA)
    
    for block_number = 1:coded_data_blocks
        if ( mod( block_number , 2 ) == 1) % odd blocks goes to 2nd Tx stream
            coded_Tx2 = [coded_Tx2  coded_data( ( N_CBPS * ( block_number - 1 ) )+1 : ( N_CBPS * block_number  ) ) ];
        else % even blocks goes to 1st Tx stream
            coded_Tx1 = [coded_Tx1 coded_data( ( N_CBPS * ( block_number - 1 ) )+1 : ( N_CBPS * block_number  ) ) ];
        end
    end
    
    
    %% Streams interleaving DATA + SIGNAL

    Modulation = T(data_rate,:).Modulation; % read modulation from T
    Modulation = Modulation{1}; % converts cell to string
    
    % interleaving 1st Tx stream
    coded_signal_length = length(coded_signal); % length of coded signal
    coded_Tx1_DATA = coded_Tx1((coded_signal_length+1):end); % coded 1st Tx stream data
    
    [idata1, idata_Tx1_DATA] = interleaver(coded_Tx1_DATA, CH_bandwidth, transmission_type, Modulation); % interleaved data 1st stream input 
    
    [idata2, idata_Tx1_SIGNAL] = interleaver(coded_signal, CH_bandwidth, transmission_type, 'BPSK'); % interleaved signal 1st stream input 
    
    % interleaving 2nd Tx stream
    
    [idata3, idata_Tx2_DATA] = interleaver(coded_Tx2, CH_bandwidth, transmission_type, Modulation); % interleaved data 2nd stream input 
    
    %% Streams subcarrier modulation DATA + SIGNAL

    % 1st stream modulation
    
    modulated_Tx1_SIGNAL = modulator(idata_Tx1_SIGNAL,'BPSK'); %modulated signal 1st stream input
    
    modulated_Tx1_DATA = modulator(idata_Tx1_DATA,Modulation); %modulated data 1st stream input
    
    % 2nd stream modulation
    
    modulated_Tx2_DATA = modulator(idata_Tx2_DATA,Modulation); %modulated data 2nd stream input
    
    %% STBC 2Txx1Rx Alamouti code
    
    h = [0.3 -.2]; %CHANNEL COEFFICENTS MATRIX
    %h11=1; h12=1; h21=1; h22=1;
    
    modulated_Tx1 = [modulated_Tx1_SIGNAL modulated_Tx1_DATA]; % 1st modulated Tx stream
    modulated_Tx2 = modulated_Tx2_DATA; % 2nd modulated Tx stream
    modulated_Tx = [modulated_Tx1 modulated_Tx2]; % modulated Tx stream 
    %power_coeff = 7.694; % power coefficient for correcting modulation
    %modulated_Tx = modulated_Tx * power_coeff;
    modulated_Tx_length = length(modulated_Tx); % length of modulated Tx stream
    Tx=zeros(modulated_Tx_length,1); 

    for i=1:(modulated_Tx_length-1)

        % Symbols at time period T;
        Tx(i,1) = modulated_Tx(i);
        Tx(i+1,1) = modulated_Tx(i+1);

        % Symbols at time period T+1;
        Tx(i,2) = -conj(modulated_Tx(i+1));
        Tx(i+1,2) = conj(modulated_Tx(i));

    end
    
    Tx1 = zeros(1,size(Tx,1));
    Tx2 = zeros(1,size(Tx,1));
    
    for i = 1:size(Tx,1)
        
        Tx1(i) = Tx(i,1); % 1st Tx stream
        Tx2(i) = Tx(i,2); % 2nd Tx stream
        
    end
    
%     Tx1 = Tx(:,1)'; % 1st Tx stream
%     Tx2 = Tx(:,2)'; % 2nd Tx stream
    

    %% Streams OFDM modulation (pitot insertion + adding GI) DATA + SIGNAL
    
    % 1st stream OFDM modulation
    [ofdm_Tx1, ifft_output_Tx1] = ofdm_modulator(Tx1,GI,1); % OFDM modulation
    
    % 2nd stream OFDM modulation
    [ofdm_Tx2, ifft_output_Tx2] = ofdm_modulator(Tx2,GI,1); % OFDM modulation
    
    ofdm_output = [ofdm_Tx1 ofdm_Tx2]; % complete DATA + SIGNAL OFDM symbols

    %% Signal power measurement before upsampling

    lFilt = 10*log10(0.8);   
    
    signal_power = 10*log10((sum(ofdm_output.*conj(ofdm_output)))/length(ofdm_output));
    signal_power = signal_power + 10*log(N_ST/BCU_W) + 2*lFilt;   % ratio of used carriers vs all carriers
    
    %% Streams pulse shaping (RRC filter) + IQ modulation
    
    % 1st stream RRC filtration + IQ modulation
    
    [I_filt_Tx1, Q_filt_Tx1] = RRC_filter(ofdm_Tx1,upsampling,beta,span); % RRC filtering
  
    [output_signal_Tx1, carrier_real_Tx1, carrier_imag_Tx1 ] = IQ_modulator(I_filt_Tx1, Q_filt_Tx1,upsampling,fc); % IQ modulation
    
    % 2nd stream RRC filtration + IQ modulation
    
    [I_filt_Tx2, Q_filt_Tx2] = RRC_filter(ofdm_Tx2,upsampling,beta,span); % RRC filtering
  
    [output_signal_Tx2, carrier_real_Tx2, carrier_imag_Tx2 ] = IQ_modulator(I_filt_Tx2, Q_filt_Tx2,upsampling,fc); % IQ modulation
    
    %% Gain
    
    power_coeff = 7.694; % power coefficient for correcting modulation
    output_signal_Tx1 = output_signal_Tx1 * power_coeff;
    output_signal_Tx2 = output_signal_Tx2 * power_coeff;
    
    perc = 40;
    waitbar(perc/100,w,sprintf('%d%% Transmitter processing...',perc));
    perc_step = floor(50/length(SNR));
    

%% RF fading channel models
   
if (result == 1) % SNR static
    
    RF_signal_Tx1 = RF_channel(selected_RF_model,output_signal_Tx1,SNR,signal_power,fs);    
    
    RF_signal_Tx2 = RF_channel(selected_RF_model,output_signal_Tx2,SNR,signal_power,fs); 
    
%     RF_signal_Tx1 = output_signal_Tx1;
%     RF_signal_Tx2 = output_signal_Tx2;
    
    RF_signal = [RF_signal_Tx1 RF_signal_Tx2];

%% Streams IQ demodulator + Filter (RRC filter)

    downsampling = upsampling; % downsampling coefficient
    
    % 1st stream RRC filtration + IQ demodulation

    [I_Rx1,Q_Rx1] = IQ_demodulator(RF_signal_Tx1, carrier_real_Tx1, carrier_imag_Tx1); % IQ demodulation

    signal_Rx1 = RRC_filter2(I_Rx1,Q_Rx1,downsampling,beta,span); % RRC filtration
    
    % 2nd stream RRC filtration + IQ demodulation

    [I_Rx2,Q_Rx2] = IQ_demodulator(RF_signal_Tx2, carrier_real_Tx2, carrier_imag_Tx2); % IQ demodulation

    signal_Rx2 = RRC_filter2(I_Rx2,Q_Rx2,downsampling,beta,span); % RRC filtration
    

%% Streams OFDM demodulator DATA + SIGNAL

    % 1st stream OFDM demodulation
    ofdm_demodulated_Rx1 = ofdm_demodulator(signal_Rx1,GI); % OFDM demodulation   
    
    % 2nd stream OFDM demodulation
    ofdm_demodulated_Rx2 = ofdm_demodulator(signal_Rx2,GI); % OFDM demodulation   
    
    ofdm_demodulated_output = [ofdm_demodulated_Rx1 ofdm_demodulated_Rx2];

%% STBC 2Txx1Rx Alamouti code

    Rx = zeros(length(ofdm_demodulated_Rx2),2); % Rx streams
    
    for i = 1:length(ofdm_demodulated_Rx1)
        
        Rx(i,1) = ofdm_demodulated_Rx1(i); % 1st Tx stream
        Rx(i,2) = ofdm_demodulated_Rx2(i); % 2nd Tx stream
        
    end
    
    modulated_Rx_length = size(Rx,1); % length of modulated Rx stream
    
    modulated_Rx = []; % modulated Rx stream
    
%     nTx = 2; % number of transmitted antennas
%     nRx = 1; % number of receiving antennas
%     hCof = zeros(2,2,modulated_Rx_length/nTx); % H^H*H [nTx * nTx] 
% 
%     hCof(1,1,:) = sum(h(:,2,:).*conj(h(:,2,:)),1);  % d term
%     hCof(2,2,:) = sum(h(:,1,:).*conj(h(:,1,:)),1);  % a term
%     hCof(2,1,:) = -sum(h(:,2,:).*conj(h(:,1,:)),1); % c term
%     hCof(1,2,:) = -sum(h(:,1,:).*conj(h(:,2,:)),1); % b term
%     
%     hDen = ((hCof(1,1,:).*hCof(2,2,:)) - (hCof(1,2,:).*hCof(2,1,:))); % ad-bc term
%     hDen = reshape(kron(reshape(hDen,1,modulated_Rx_length/nTx),ones(2,2)),2,2,modulated_Rx_length/nTx);  % formatting for division
%     hInv = hCof./hDen; % inv(H^H*H)
%     
%     hMod =  reshape(conj(h),nRx,modulated_Rx_length); % H^H operation
%     
%     yMod = kron(ofdm_demodulated_Rx2,ones(1,2)); % formatting the received symbol for equalization
%     yMod = sum(hMod.*yMod,1); % H^H * y 
%     yMod =  kron(reshape(yMod,2,modulated_Rx_length/nTx),ones(1,2)); % formatting
%     yHat = sum(reshape(hInv,2,modulated_Rx_length).*yMod,1); % inv(H^H*H)*H^H*y
%     
%     % receiver - hard decision decoding
%     ipHat = real(yHat)>0;
    
    for i=1:2:modulated_Rx_length
        
        s1=Rx(i,1);
        s2=Rx(i+1,1);

        %Recieved data by RX1 Antenna at time interval T
        r(1,1)= (h(1,1)*s1) + (h(1,2)*s2); % + e(1,1);

        %Recieved data by RX1 Antenna at time interval (T+1)
        r(1,2)= ((-h(1,1))*conj(s2)) + (h(1,2)*conj(s1)); % + e(1,2);

        t(1,1)=((conj(h(1,1))*r(1,1)));
        t(1,2)=h(1,2)*(conj(r(1,2)));
        t(2,1)=((conj(h(1,2)))*r(1,1));
        t(2,2)=((h(1,1)*(conj(r(1,2)))));

        %Maximum Likelehhod Detection Scehme
        s1_e =t(1,1) + t(1,2);
        s2_e= t(2,1) - t(2,2);

        Rx_per_period = [s1_e, s2_e];

        modulated_Rx = [modulated_Rx Rx_per_period];
    
    end
    % power_coeff = 7.694; % power coefficient for correcting modulation
    % modulated_Rx = modulated_Rx * power_coeff; 
    
    N_CBPS = 48; % number of coded bits per symbol in BPSK modulation
    modulated_Rx1_SIGNAL = modulated_Rx(1:N_CBPS); % 1st modulated Rx stream (SIGNAL)
    modulated_Rx_DATA = modulated_Rx(N_CBPS+1:end); % Rx stream DATA
    modulated_Rx_DATA_length = length(modulated_Rx_DATA); % length of Rx stream DATA
    
    N_BPSC = T(data_rate,:).N_BPSC; % read number of coded bits per subcarrier from T
    N_CBPS = T(data_rate,:).N_CBPS; % read number of coded bits per symbol from T
    coded_DATA_Rx_length = modulated_Rx_DATA_length * N_BPSC; % number of coded Rx DATA bits
    number_of_blocks = coded_DATA_Rx_length / N_CBPS; % number of blocks
    
    % Recognition of DATA blocks contribution between streams
    
    coded_Rx1 = modulated_Rx1_SIGNAL; % 1st coded Rx stream (SIGNAL + DATA) 
    coded_Rx2 = []; % 2nd coded Rx stream (DATA)
    
    coded_data = zeros(1,coded_DATA_Rx_length);
    
    for block_number = 1:number_of_blocks
        if ( mod( block_number , 2 ) == 1) % odd blocks goes to 2nd Tx stream
            coded_Rx2 = [coded_Rx2  coded_data( ( N_CBPS * ( block_number - 1 ) )+1 : ( N_CBPS * block_number  ) ) ];
        else % even blocks goes to 1st Tx stream
            coded_Rx1 = [coded_Rx1 coded_data( ( N_CBPS * ( block_number - 1 ) )+1 : ( N_CBPS * block_number  ) ) ];
        end
    end
    
    modulated_Rx1_SIGNAL_length = length(modulated_Rx1_SIGNAL); % length of SIGNAL
    modulated_Rx1_DATA_length = ( length(coded_Rx1) - modulated_Rx1_SIGNAL_length ) / N_BPSC;
    modulated_Rx2_DATA_length = length(coded_Rx2) / N_BPSC;

    modulated_Rx1_DATA = modulated_Rx_DATA( 1 : modulated_Rx1_DATA_length ); % 1st Rx stream (DATA)
    modulated_Rx2_DATA = modulated_Rx_DATA( modulated_Rx1_DATA_length + 1 : end ); % 2nd Rx stream (DATA)
    

%% Streams subcarrier demodulator DATA + SIGNAL
    
    % 1st stream demodulation
    
    demodulated_Rx1_SIGNAL = demodulator(modulated_Rx1_SIGNAL,'BPSK'); % demodulated 1st Rx stream (SIGNAL)   
    demodulated_Rx1_DATA = demodulator(modulated_Rx1_DATA,Modulation); % demodulated 1st Rx stream (DATA)
    
    % 2nd stream demodulation
    
    demodulated_Rx2_DATA = demodulator(modulated_Rx2_DATA,Modulation); % demodulated 2nd Rx stream (DATA)
    
    

%% Streams deinterleaving DATA + SIGNAL
    
    % deinterleaving 1st Rx stream
    
    deinterleaved_Rx1_SIGNAL = deinterleaver(demodulated_Rx1_SIGNAL, CH_bandwidth, transmission_type, 'BPSK'); % deinterleaved 1st Rx stream (SIGNAL)
    deinterleaved_Rx1_DATA = deinterleaver(demodulated_Rx1_DATA, CH_bandwidth, transmission_type, Modulation); % deinterleaved 1st Rx stream (DATA) 
    
    % deinterleaving 2nd Rx stream
    
    deinterleaved_Rx2_DATA = deinterleaver(demodulated_Rx2_DATA, CH_bandwidth, transmission_type, Modulation); % deinterleaved 2nd Rx stream (DATA) 
    
    deinterleaved_signal = deinterleaved_Rx1_SIGNAL;
    deinterleaved_data = [deinterleaved_Rx1_DATA deinterleaved_Rx2_DATA];
    
%% FEC decoding + puncturing DATA + SIGNAL

    decoded_signal = fec_decoder(deinterleaved_signal,'1/2'); % convolutionally decoded signal

    decoded_data = fec_decoder(deinterleaved_data,R); % convolutionally decoded data

%% Descrambler DATA
    
    descramblered_data = descrambler(decoded_data,s,half_byte); % descrambler DATA

%% BER calculation
    
    TxData = [coded_signal coded_data]; % transmitted bits
    
    % pre Viterbi decoder
    RxData = [deinterleaved_signal deinterleaved_data]; % received bits before Viterbi decoder
    nErrors = biterr(TxData,RxData); % calculate the number of bit errors
    numErrs = nErrors; % number of error bits
    numBits = length(TxData); % number of transmitted bits
    pre_BER_vector = numErrs / numBits; % BER calculated
    
    TxData = [SIGNAL scrambled_data]; % transmitted bits
    
    % post Viterbi decoder
    RxData = [decoded_signal descramblered_data]; % received bits after Viterbi decoder
    nErrors = biterr(TxData,RxData); % calculate the number of bit errors
    numErrs = nErrors; % number of error bits
    numBits = length(TxData); % number of transmitted bits
    post_BER_vector = numErrs / numBits; % BER calculated
    
    perc = perc + perc_step;
    waitbar(perc/100,w,sprintf('%d%% Receiver processing...',perc));
    
%% MER + EVM calculation

    MER_vector = MER_calc(modulated_Tx, ofdm_demodulated_output); % MER calculation
    
    EVM_vector = EVM_calc(modulated_Tx, ofdm_demodulated_output); % EVM calculation
    
    modulated_output = modulated_Tx;
    
else % SNR dynamic
    
    pre_BER_vector = []; % vector for BER before Viterbi decoder on each SNR
    post_BER_vector = []; % vector for BER after Viterbi decoder on each SNR
    MER_vector = []; % vector for MER calculations on each SNR
    EVM_vector = []; % vector for EVM calculations on each SNR
    pre_TxData = [coded_signal coded_data]; % transmitted bits
    post_TxData = [SIGNAL scrambled_data]; % transmitted bits
    pre_numBits = length(pre_TxData); % number of transmitted bits ( 1x transmission )
    post_numBits = length(post_TxData); % number of transmitted bits ( 1x transmission )  
    
    for n=1:length(SNR)
        
        RxData = 0;
        
        RF_signal_Tx1 = RF_channel(selected_RF_model,output_signal_Tx1,SNR(n),signal_power,fs);    
    
        RF_signal_Tx2 = RF_channel(selected_RF_model,output_signal_Tx2,SNR(n),signal_power,fs); 

        RF_signal = [RF_signal_Tx1 RF_signal_Tx2];
        
        %% IQ demodulator + Filter (RRC filter)

        downsampling = upsampling; % downsampling coefficient

        % 1st stream RRC filtration + IQ demodulation

        [I_Rx1,Q_Rx1] = IQ_demodulator(RF_signal_Tx1, carrier_real_Tx1, carrier_imag_Tx1); % IQ demodulation

         signal_Rx1 = RRC_filter2(I_Rx1,Q_Rx1,downsampling,beta,span); % RRC filtration
    
        % 2nd stream RRC filtration + IQ demodulation

        [I_Rx2,Q_Rx2] = IQ_demodulator(RF_signal_Tx2, carrier_real_Tx2, carrier_imag_Tx2); % IQ demodulation

        signal_Rx2 = RRC_filter2(I_Rx2,Q_Rx2,downsampling,beta,span); % RRC filtration
    
        %% Streams OFDM demodulator DATA + SIGNAL

        % 1st stream OFDM demodulation
        ofdm_demodulated_Rx1 = ofdm_demodulator(signal_Rx1,GI); % OFDM demodulation   

        % 2nd stream OFDM demodulation
        ofdm_demodulated_Rx2 = ofdm_demodulator(signal_Rx2,GI); % OFDM demodulation   

        ofdm_demodulated_output = [ofdm_demodulated_Rx1 ofdm_demodulated_Rx2];

        %% STBC 2Txx1Rx Alamouti code

        Rx = zeros(length(ofdm_demodulated_Rx2),2); % Rx streams

        for i = 1:length(ofdm_demodulated_Rx1)

            Rx(i,1) = ofdm_demodulated_Rx1(i); % 1st Tx stream
            Rx(i,2) = ofdm_demodulated_Rx2(i); % 2nd Tx stream

        end

        modulated_Rx_length = size(Rx,1); % length of modulated Rx stream

        modulated_Rx = []; % modulated Rx stream

        for i=1:2:modulated_Rx_length

            s1=Rx(i,1);
            s2=Rx(i+1,1);

            %Recieved data by RX1 Antenna at time interval T
            r(1,1)= (h(1,1)*s1) + (h(1,2)*s2); % + e(1,1);

            %Recieved data by RX1 Antenna at time interval (T+1)
            r(1,2)= ((-h(1,1))*conj(s2)) + (h(1,2)*conj(s1)); % + e(1,2);

            t(1,1)=((conj(h(1,1))*r(1,1)));
            t(1,2)=h(1,2)*(conj(r(1,2)));
            t(2,1)=((conj(h(1,2)))*r(1,1));
            t(2,2)=((h(1,1)*(conj(r(1,2)))));

            %Maximum Likelehhod Detection Scehme
            s1_e =t(1,1) + t(1,2);
            s2_e= t(2,1) - t(2,2);

            Rx_per_period = [s1_e, s2_e];

            modulated_Rx = [modulated_Rx Rx_per_period];

        end
        % power_coeff = 7.694; % power coefficient for correcting modulation
        % modulated_Rx = modulated_Rx * power_coeff; 

        N_CBPS = 48; % number of coded bits per symbol in BPSK modulation
        modulated_Rx1_SIGNAL = modulated_Rx(1:N_CBPS); % 1st modulated Rx stream (SIGNAL)
        modulated_Rx_DATA = modulated_Rx(N_CBPS+1:end); % Rx stream DATA
        modulated_Rx_DATA_length = length(modulated_Rx_DATA); % length of Rx stream DATA

        N_BPSC = T(data_rate,:).N_BPSC; % read number of coded bits per subcarrier from T
        N_CBPS = T(data_rate,:).N_CBPS; % read number of coded bits per symbol from T
        coded_DATA_Rx_length = modulated_Rx_DATA_length * N_BPSC; % number of coded Rx DATA bits
        number_of_blocks = coded_DATA_Rx_length / N_CBPS; % number of blocks

        % Recognition of DATA blocks contribution between streams

        coded_Rx1 = modulated_Rx1_SIGNAL; % 1st coded Rx stream (SIGNAL + DATA) 
        coded_Rx2 = []; % 2nd coded Rx stream (DATA)

        coded_data = zeros(1,coded_DATA_Rx_length);

        for block_number = 1:number_of_blocks
            if ( mod( block_number , 2 ) == 1) % odd blocks goes to 2nd Tx stream
                coded_Rx2 = [coded_Rx2  coded_data( ( N_CBPS * ( block_number - 1 ) )+1 : ( N_CBPS * block_number  ) ) ];
            else % even blocks goes to 1st Tx stream
                coded_Rx1 = [coded_Rx1 coded_data( ( N_CBPS * ( block_number - 1 ) )+1 : ( N_CBPS * block_number  ) ) ];
            end
        end

        modulated_Rx1_SIGNAL_length = length(modulated_Rx1_SIGNAL); % length of SIGNAL
        modulated_Rx1_DATA_length = ( length(coded_Rx1) - modulated_Rx1_SIGNAL_length ) / N_BPSC;
        modulated_Rx2_DATA_length = length(coded_Rx2) / N_BPSC;

        modulated_Rx1_DATA = modulated_Rx_DATA( 1 : modulated_Rx1_DATA_length ); % 1st Rx stream (DATA)
        modulated_Rx2_DATA = modulated_Rx_DATA( modulated_Rx1_DATA_length + 1 : end ); % 2nd Rx stream (DATA)


        %% Streams subcarrier demodulator DATA + SIGNAL

        % 1st stream demodulation

        demodulated_Rx1_SIGNAL = demodulator(modulated_Rx1_SIGNAL,'BPSK'); % demodulated 1st Rx stream (SIGNAL)   
        demodulated_Rx1_DATA = demodulator(modulated_Rx1_DATA,Modulation); % demodulated 1st Rx stream (DATA)

        % 2nd stream demodulation

        demodulated_Rx2_DATA = demodulator(modulated_Rx2_DATA,Modulation); % demodulated 2nd Rx stream (DATA)



        %% Streams deinterleaving DATA + SIGNAL

        % deinterleaving 1st Rx stream

        deinterleaved_Rx1_SIGNAL = deinterleaver(demodulated_Rx1_SIGNAL, CH_bandwidth, transmission_type, 'BPSK'); % deinterleaved 1st Rx stream (SIGNAL)
        deinterleaved_Rx1_DATA = deinterleaver(demodulated_Rx1_DATA, CH_bandwidth, transmission_type, Modulation); % deinterleaved 1st Rx stream (DATA) 

        % deinterleaving 2nd Rx stream

        deinterleaved_Rx2_DATA = deinterleaver(demodulated_Rx2_DATA, CH_bandwidth, transmission_type, Modulation); % deinterleaved 2nd Rx stream (DATA) 

        deinterleaved_signal = deinterleaved_Rx1_SIGNAL;
        deinterleaved_data = [deinterleaved_Rx1_DATA deinterleaved_Rx2_DATA];

        %% FEC decoding + puncturing DATA + SIGNAL

        decoded_signal = fec_decoder(deinterleaved_signal,'1/2'); % convolutionally decoded signal

        decoded_data = fec_decoder(deinterleaved_data,R); % convolutionally decoded data

        %% Descrambler DATA

        descramblered_data = descrambler(decoded_data,s,half_byte); % descrambler DATA
        
        %% BER calculation
        
        TxData = [coded_signal coded_data]; % transmitted bits
        
        % pre Viterbi decoder
        RxData = [deinterleaved_signal deinterleaved_data]; % received bits before Viterbi decoder
        nErrors = biterr(pre_TxData,RxData); % calculate the number of bit errors
        numErrs = nErrors; % number of error bits
        %numBits = length(pre_TxData); % number of transmitted bits
        pre_BER_actual = numErrs / pre_numBits; % BER before Viterbi decoder calculated
        pre_BER_vector = [pre_BER_vector pre_BER_actual]; % fill BER vector with calculated BER before Viterbi decoder on each SNR
    
        % post Viterbi decoder
        RxData = [decoded_signal decoded_data]; % received bits after Viterbi decoder
        nErrors = biterr(post_TxData,RxData); % calculate the number of bit errors
        numErrs = nErrors; % number of error bits
        %numBits = length(TxData); % number of transmitted bits
        post_BER_actual = numErrs / post_numBits; % BER after Viterbi decoder calculated
        post_BER_vector = [post_BER_vector post_BER_actual]; % fill BER vector with calculated BER after Viterbi decoder on each SNR
        
        perc = perc + perc_step;
        waitbar(perc/100,w,sprintf('%d%% Receiver processing...',perc));
        
        %% MER + EVM calculation

        MER_act = MER_calc(modulated_Tx, ofdm_demodulated_output); % MER calculation
        MER_vector = [MER_vector MER_act]; % MER vector accumulation

        EVM_act = EVM_calc(modulated_Tx, ofdm_demodulated_output); % EVM calculation
        EVM_vector = [EVM_vector EVM_act]; % EVM vector accumulation
        
    end
    modulated_output = modulated_Tx;
end

elseif ( ( N_Tx == 2 ) && ( N_Rx == 2 ) ) % if used MIMO mode
    
    %% Stream parser 
    
    N_CBPS = T(MCS_index,:).N_CBPS; % read number of coded bits per symbol from T
    coded_data_blocks = length(coded_data) / N_CBPS; % number of coded data blocks
    
    %coded_signal for non-HT
    %coded_data for VHT - does not work
    
    %coded_signal = fec_coder(SIGNAL,'1/2'); % convolutionally coded signal input  
    coded_Tx1 = coded_signal; % coded_signal; % 1st coded Tx stream (SIGNAL + DATA) 
    coded_Tx2 = []; % 2nd coded Tx stream (DATA)
    
    for block_number = 1:coded_data_blocks
        if ( mod( block_number , 2 ) == 1) % odd blocks goes to 2nd Tx stream
            coded_Tx2 = [coded_Tx2  coded_data( ( N_CBPS * ( block_number - 1 ) )+1 : ( N_CBPS * block_number  ) ) ];
        else % even blocks to 1st Tx stream
            coded_Tx1 = [coded_Tx1 coded_data( ( N_CBPS * ( block_number - 1 ) )+1 : ( N_CBPS * block_number  ) ) ];
        end
    end
    
    
    %% Streams interleaving DATA + SIGNAL

    Modulation = T(MCS_index,:).Modulation; % read modulation from T
    Modulation = Modulation{1}; % converts cell to string
    
    % interleaving 1st Tx stream
    coded_signal_length = length(coded_signal); %coded_signal; % length of coded signal
    coded_Tx1_DATA = coded_Tx1((coded_signal_length+1):end); % coded 1st Tx stream data
    
    [idata1, idata_Tx1_DATA] = interleaver(coded_Tx1_DATA, CH_bandwidth, transmission_type, Modulation); % interleaved data 1st stream input 
    
    [idata2, idata_Tx1_SIGNAL] = interleaver(coded_signal, CH_bandwidth, transmission_type,'BPSK'); %coded_signal; % interleaved signal 1st stream input 
    
    % interleaving 2nd Tx stream
    
    [idata3, idata_Tx2_DATA] = interleaver(coded_Tx2, CH_bandwidth, transmission_type,Modulation); % interleaved data 2nd stream input 
    
    %% Streams subcarrier modulation DATA + SIGNAL

    % 1st stream modulation
    
    modulated_Tx1_SIGNAL = modulator(idata_Tx1_SIGNAL,'BPSK'); %modulated signal 1st stream input
    
    modulated_Tx1_DATA = modulator(idata_Tx1_DATA,Modulation); %modulated data 1st stream input
    
    % 2nd stream modulation
    
    modulated_Tx2_DATA = modulator(idata_Tx2_DATA,Modulation); %modulated data 2nd stream input
    
    %% STBC 2Txx2Rx Alamouti code
    
    h=[.5 .5 ;.5 .5]; %CHANNEL COEFFICENTS MATRIX
    %h=[1 1 ;1 1]; %CHANNEL COEFFICENTS MATRIX
    %h11=1; h12=1; h21=1; h22=1;
    
    modulated_Tx1 = [modulated_Tx1_SIGNAL modulated_Tx1_DATA]; % 1st modulated Tx stream
    modulated_Tx2 = modulated_Tx2_DATA; % 2nd modulated Tx stream
    modulated_Tx = [modulated_Tx1 modulated_Tx2]; % modulated Tx stream 
    %power_coeff = 7.694; % power coefficient for correcting modulation
    %modulated_Tx = modulated_Tx * power_coeff;
    modulated_Tx_length = length(modulated_Tx); % length of modulated Tx stream
    Tx=zeros(modulated_Tx_length,1); 

    for i=1:(modulated_Tx_length-1)

        % Symbols at time period T;
        Tx(i,1) = h(1,1) * modulated_Tx(i);
        Tx(i+1,1) = h(2,1) * modulated_Tx(i+1);

        % Symbols at time period T+1;
        Tx(i,2) = -conj(modulated_Tx(i+1));
        Tx(i+1,2) = conj(modulated_Tx(i));

    end
    
    Tx1 = zeros(1,size(Tx,1));
    Tx2 = zeros(1,size(Tx,1));
    
    for i = 1:size(Tx,1)
        
        Tx1(i) = Tx(i,1); % 1st Tx stream
        Tx2(i) = Tx(i,2); % 2nd Tx stream
        
    end
    
%    Tx1 = Tx(:,1)'; % 1st Tx stream
%    Tx2 = Tx(:,2)'; % 2nd Tx stream
    

    %% Streams OFDM modulation (pitot insertion + adding GI) DATA + SIGNAL
    
    % 1st stream OFDM modulation
    [ofdm_Tx1, Tx1_flag] = ofdm_modulator2(Tx1,transmission_type,CH_bandwidth,BCU_bandwidth,GI); % OFDM modulation

    % 2nd stream OFDM modulation
    [ofdm_Tx2, Tx2_flag] = ofdm_modulator2(Tx2,transmission_type,CH_bandwidth,BCU_bandwidth,GI); % OFDM modulation

    ofdm_output = [ofdm_Tx1 ofdm_Tx2]; % complete DATA + SIGNAL OFDM symbols

    %% Signal power measurement before upsampling

    %lFilt = 10*log10(0.8);   
    
    %signal_power = 10*log10((sum(ofdm_output.*conj(ofdm_output)))/length(ofdm_output));
    %signal_power = signal_power + 10*log(N_ST/BCU_W) + 2*lFilt;   % ratio of used carriers vs all carriers
    
    lFilt = 10*log10(0.8);   
    
    signal_power = 10*log10((sum(ofdm_output.*conj(ofdm_output)))/length(ofdm_output));
    %signal_power = signal_power + 10*log(N_ST/BCU_W) + 2*lFilt;   % ratio of used carriers vs all carriers
    signal_power = signal_power + 10*log((N_SD-N_SP)/BCU_W) + 2*lFilt;   % ratio of used carriers vs all carriers
    
    %% Streams pulse shaping (RRC filter) + IQ modulation
    
    % 1st stream RRC filtration + IQ modulation
    
    [I_filt_Tx1, Q_filt_Tx1] = RRC_filter(ofdm_Tx1,upsampling,beta,span); % RRC filtering
  
    [output_signal_Tx1, carrier_real_Tx1, carrier_imag_Tx1 ] = IQ_modulator(I_filt_Tx1, Q_filt_Tx1,upsampling,fc, delta_f); % IQ modulation
    %[I_Rx1,Q_Rx1] = IQ_demodulator(output_signal_Tx1, carrier_real_Tx1, carrier_imag_Tx1); % IQ demodulation
    % 2nd stream RRC filtration + IQ modulation
    
    [I_filt_Tx2, Q_filt_Tx2] = RRC_filter(ofdm_Tx2,upsampling,beta,span); % RRC filtering
  
    [output_signal_Tx2, carrier_real_Tx2, carrier_imag_Tx2 ] = IQ_modulator(I_filt_Tx2, Q_filt_Tx2,upsampling,fc, delta_f); % IQ modulation
    
    %% Gain
    
    %power_coeff = 7.694; % power coefficient for correcting modulation
    %output_signal_Tx1 = output_signal_Tx1 * power_coeff;
    %output_signal_Tx2 = output_signal_Tx2 * power_coeff;
    
     perc = 40;
     waitbar(perc/100,w,sprintf('%d%% Transmitter processing...',perc));
     perc_step = floor(50/length(SNR));
    

%% RF fading channel models
   
if (result == 1) % SNR static
    
    RF_signal_Tx1 = RF_channel(selected_RF_model,output_signal_Tx1,SNR,signal_power,fs);    
    
    RF_signal_Tx2 = RF_channel(selected_RF_model,output_signal_Tx2,SNR,signal_power,fs); 
    
    %RF_signal_Tx1 = output_signal_Tx1;
    %RF_signal_Tx2 = output_signal_Tx2;
    
    RF_signal = [RF_signal_Tx1 RF_signal_Tx2];

%% Streams IQ demodulator + Filter (RRC filter)

    downsampling = upsampling; % downsampling coefficient
    
    % 1st stream RRC filtration + IQ demodulation

    [I_Rx1,Q_Rx1] = IQ_demodulator(RF_signal_Tx1, carrier_real_Tx1, carrier_imag_Tx1); % IQ demodulation

    signal_Rx1 = RRC_filter2(I_Rx1,Q_Rx1,downsampling,beta,span); % RRC filtration
    
    % 2nd stream RRC filtration + IQ demodulation

    [I_Rx2,Q_Rx2] = IQ_demodulator(RF_signal_Tx2, carrier_real_Tx2, carrier_imag_Tx2); % IQ demodulation

    signal_Rx2 = RRC_filter2(I_Rx2,Q_Rx2,downsampling,beta,span); % RRC filtration
    

%% Streams OFDM demodulator DATA + SIGNAL

    % 1st stream OFDM demodulation
    [ofdm_demodulated_Rx1, ofdm_demodulated_Rx1_pilots] = ofdm_demodulator2(signal_Rx1,transmission_type,CH_bandwidth,BCU_bandwidth,GI,Tx1_flag,equalizer_flag); % OFDM demodulation
    
    % 2nd stream OFDM demodulation
    [ofdm_demodulated_Rx2, ofdm_demodulated_Rx2_pilots] = ofdm_demodulator2(signal_Rx2,transmission_type,CH_bandwidth,BCU_bandwidth,GI,Tx2_flag,equalizer_flag); % OFDM demodulation
    
    ofdm_demodulated_output = [ofdm_demodulated_Rx1 ofdm_demodulated_Rx2];
    
    ofdm_demodulated_output_pilots = [ofdm_demodulated_Rx1_pilots ofdm_demodulated_Rx2_pilots];

%% STBC 2Txx2Rx Alamouti code
%ofdm_demodulated_Rx1 = Tx1;
%ofdm_demodulated_Rx2 = Tx2;
    Rx = zeros(length(ofdm_demodulated_Rx2),2); % Rx streams
    
    for i = 1:length(ofdm_demodulated_Rx1)     
        Rx(i,1) = ofdm_demodulated_Rx1(i); % 1st Tx stream
        Rx(i,2) = ofdm_demodulated_Rx2(i); % 2nd Tx stream
%         Rx(i,1) = Tx1(i); % 1st Tx stream
%         Rx(i,2) = Tx2(i); % 2nd Tx stream        
    end
    
    modulated_Rx_length = size(Rx,1); % length of modulated Rx stream
    
    modulated_Rx = []; % modulated Rx stream
    
    for i=1:2:modulated_Rx_length
        
        s1=Rx(i,1);
        s2=Rx(i+1,1);

        %Recieved data by RX1 Antenna at time interval T
        r(1,1)= (h(1,1)*s1) + (h(1,2)*s2); % + e(1,1);

        %Recieved data by RX1 Antenna at time interval (T+1)
        r(1,2)= ((-h(1,1))*conj(s2)) + (h(1,2)*conj(s1)); % + e(1,2);

        %Recieved data by RX2 Antenna at time interval T
        r(2,1)= (h(2,1)*s1) + (h(2,2)*s2); % + e(2,1);

        %Recieved data by RX2 Antenna at time interval (T+1)
        r(2,2)= ((-h(2,1))*conj(s2)) + (h(2,2)*conj(s1)); % + e(2,2);
        
        t(1,1)=((conj(h(1,1))*r(1,1)));
        t(1,2)=h(1,2)*(conj(r(1,2)));
        t(2,1)=((conj(h(2,1)))*r(2,1));
        t(2,2)=((h(1,2)*(conj(r(2,2)))));
        
        % channel matrix aproximation
        c(1,1)= ((conj(h(1,2)))*r(1,1));
        c(1,2)= h(1,1)*(conj(r(1,2)));
        c(2,1)= ((conj(h(2,2)))*r(2,1));
        c(2,2)= ((h(2,1)*(conj(r(2,2)))));

        %Maximum Likelehhod Detection Scehme
        s1_e =t(1,1) + t(1,2) + t(2,1) + t(2,2);
        s2_e= c(1,1) - c(1,2) + c(2,1) - c(2,2);

        Rx_per_period = [s1_e, s2_e];

        modulated_Rx = [modulated_Rx Rx_per_period];
    
    end
    % power_coeff = 7.694; % power coefficient for correcting modulation
    % modulated_Rx = modulated_Rx * power_coeff; 
    
%     N_CBPS_signal = 48; % number of coded bits per symbol in BPSK modulation
%     
%     modulated_Tx1 = [modulated_Tx1_SIGNAL modulated_Tx1_DATA]; % 1st modulated Tx stream
%     modulated_Tx2 = modulated_Tx2_DATA; % 2nd modulated Tx stream
%     modulated_Tx = [modulated_Tx1 modulated_Tx2]; % modulated Tx stream 
%     %power_coeff = 7.694; % power coefficient for correcting modulation
%     %modulated_Tx = modulated_Tx * power_coeff;
%     modulated_Tx_length = length(modulated_Tx); % length of modulated Tx stream

    N_CBPS = 48; % number of coded bits per symbol in BPSK modulation
    modulated_Rx1_SIGNAL = modulated_Rx(1:N_CBPS); % 1st modulated Rx stream (SIGNAL)
    modulated_Rx_DATA = modulated_Rx(N_CBPS+1:end); % Rx stream DATA
    modulated_Rx_DATA_length = length(modulated_Rx_DATA); % length of Rx stream DATA
    
    N_BPSC = T(MCS_index,:).N_BPSC; % read number of coded bits per subcarrier from T
    N_CBPS = T(MCS_index,:).N_CBPS; % read number of coded bits per symbol from T
    coded_DATA_Rx_length = modulated_Rx_DATA_length * N_BPSC; % number of coded Rx DATA bits
    number_of_blocks = coded_DATA_Rx_length / N_CBPS; % number of blocks
    
    % Recognition of DATA blocks contribution between streams
    
    coded_Rx1 = modulated_Rx1_SIGNAL; % 1st coded Rx stream (SIGNAL + DATA) 
    coded_Rx2 = []; % 2nd coded Rx stream (DATA)
    
    coded_data_zeros = zeros(1,coded_DATA_Rx_length);
    
    for block_number = 1:number_of_blocks
        if ( mod( block_number , 2 ) == 1) % odd blocks goes to 2nd Tx stream
            coded_Rx2 = [coded_Rx2  coded_data_zeros( ( N_CBPS * ( block_number - 1 ) )+1 : ( N_CBPS * block_number  ) ) ];
        else % even blocks goes to 1st Tx stream
            coded_Rx1 = [coded_Rx1 coded_data_zeros( ( N_CBPS * ( block_number - 1 ) )+1 : ( N_CBPS * block_number  ) ) ];
        end
    end
    
    modulated_Rx1_SIGNAL_length = length(modulated_Rx1_SIGNAL); % length of SIGNAL
    modulated_Rx1_DATA_length = ( length(coded_Rx1) - modulated_Rx1_SIGNAL_length ) / N_BPSC;
    modulated_Rx2_DATA_length = length(coded_Rx2) / N_BPSC;

    modulated_Rx1_DATA = modulated_Rx_DATA( 1 : modulated_Rx1_DATA_length ); % 1st Rx stream (DATA)
    modulated_Rx2_DATA = modulated_Rx_DATA( modulated_Rx1_DATA_length + 1 : end ); % 2nd Rx stream (DATA)
    

%% Streams subcarrier demodulator DATA + SIGNALd
    
    % 1st stream demodulation
    
    demodulated_Rx1_SIGNAL = demodulator(modulated_Rx1_SIGNAL,'BPSK'); % demodulated 1st Rx stream (SIGNAL)   
    demodulated_Rx1_DATA = demodulator(modulated_Rx1_DATA,Modulation); % demodulated 1st Rx stream (DATA)
    
    % 2nd stream demodulation
    
    demodulated_Rx2_DATA = demodulator(modulated_Rx2_DATA,Modulation); % demodulated 2nd Rx stream (DATA)

%% Streams deinterleaving DATA + SIGNAL
    
    % deinterleaving 1st Rx stream
    
    deinterleaved_Rx1_SIGNAL = deinterleaver(demodulated_Rx1_SIGNAL, CH_bandwidth, transmission_type, 'BPSK'); % deinterleaved 1st Rx stream (SIGNAL)
    deinterleaved_Rx1_DATA = deinterleaver(demodulated_Rx1_DATA, CH_bandwidth, transmission_type, Modulation); % deinterleaved 1st Rx stream (DATA) 
    
    % deinterleaving 2nd Rx stream
    
    deinterleaved_Rx2_DATA = deinterleaver(demodulated_Rx2_DATA, CH_bandwidth, transmission_type, Modulation); % deinterleaved 2nd Rx stream (DATA) 
    
    %deinterleaved_signal = deinterleaved_Rx1_SIGNAL;
    %deinterleaved_data = [deinterleaved_Rx1_DATA deinterleaved_Rx2_DATA];

%% Steam parser

    N_CBPS = T(MCS_index,:).N_CBPS; % read number of coded bits per symbol from T
    coded_data_blocks_Rx1 = length(deinterleaved_Rx1_DATA) / N_CBPS; % number of coded data blocks Rx1
    coded_data_blocks_Rx2 = length(deinterleaved_Rx2_DATA) / N_CBPS; % number of coded data blocks Rx2
    coded_data_blocks = ( length(deinterleaved_Rx1_DATA) + length(deinterleaved_Rx2_DATA) ) / N_CBPS; % number of coded data blocks

    Rx1_DATA = reshape(deinterleaved_Rx1_DATA.',N_CBPS,coded_data_blocks_Rx1).';
    Rx2_DATA = reshape(deinterleaved_Rx2_DATA.',N_CBPS,coded_data_blocks_Rx2).';
  
    deinterleaved_data = []; % coded data
    j1 = 1; j2 = 1;
    
    for block_number = 1:coded_data_blocks
        if ( mod( block_number , 2 ) == 1) % odd blocks get from data stream Rx2
            deinterleaved_data = [deinterleaved_data Rx2_DATA(j1,:)];
            j1 = j1 + 1;
        else % even blocks get from data stream Rx1
            deinterleaved_data = [deinterleaved_data Rx1_DATA(j2,:)];
            j2 = j2 + 1;
        end
    end    
    
    deinterleaved_signal = deinterleaved_Rx1_SIGNAL;
    
%% FEC decoding + puncturing DATA + SIGNAL

    decoded_signal = fec_decoder(deinterleaved_signal,'1/2'); % convolutionally decoded signal

    decoded_data = fec_decoder(deinterleaved_data,R); % convolutionally decoded data

%% Descrambler DATA
    
    descramblered_data = descrambler(decoded_data,s,half_byte); % descrambler DATA

%% BER calculation
    
    TxData = [coded_signal coded_data]; % transmitted bits
    
    % pre Viterbi decoder
    RxData = [deinterleaved_signal deinterleaved_data]; % received bits before Viterbi decoder
    nErrors = biterr(TxData,RxData); % calculate the number of bit errors
    numErrs = nErrors; % number of error bits
    numBits = length(TxData); % number of transmitted bits
    pre_BER_vector = numErrs / numBits; % BER calculated
    
    TxData = [SIGNAL scrambled_data]; % transmitted bits
    %TxData = [SIGNAL DATA]; % transmitted bits
    
    % post Viterbi decoder
    RxData = [decoded_signal decoded_data]; % received bits after Viterbi decoder
    %RxData = [decoded_signal descramblered_data]; % received bits after Viterbi decoder
    nErrors = biterr(TxData,RxData); % calculate the number of bit errors
    numErrs = nErrors; % number of error bits
    numBits = length(TxData); % number of transmitted bits
    post_BER_vector = numErrs / numBits; % BER calculated
    
    perc = perc + perc_step;
    waitbar(perc/100,w,sprintf('%d%% Receiver processing...',perc));
    
%% MER + EVM calculation

    MER_vector = MER_calc(modulated_Tx, ofdm_demodulated_output); % MER calculation
    
    EVM_vector = EVM_calc(modulated_Tx, ofdm_demodulated_output); % EVM calculation
    
    modulated_output = modulated_Tx;
    
else % SNR dynamic
    
    pre_BER_vector = []; % vector for BER before Viterbi decoder on each SNR
    post_BER_vector = []; % vector for BER after Viterbi decoder on each SNR
    MER_vector = []; % vector for MER calculations on each SNR
    EVM_vector = []; % vector for EVM calculations on each SNR
    pre_TxData = [coded_signal coded_data]; % transmitted bits
    post_TxData = [SIGNAL scrambled_data]; % transmitted bits
    pre_numBits = length(pre_TxData); % number of transmitted bits ( 1x transmission )
    post_numBits = length(post_TxData); % number of transmitted bits ( 1x transmission )  
    
    for n=1:length(SNR)
        
        pre_RxData = 0;
        post_RxData = 0;        
        RxData = 0;
        
        RF_signal_Tx1 = RF_channel(selected_RF_model,output_signal_Tx1,SNR(n),signal_power,fs);    
    
        RF_signal_Tx2 = RF_channel(selected_RF_model,output_signal_Tx2,SNR(n),signal_power,fs); 

        RF_signal = [RF_signal_Tx1 RF_signal_Tx2];
        
        %% IQ demodulator + Filter (RRC filter)

        downsampling = upsampling; % downsampling coefficient

        % 1st stream RRC filtration + IQ demodulation

        [I_Rx1,Q_Rx1] = IQ_demodulator(RF_signal_Tx1, carrier_real_Tx1, carrier_imag_Tx1); % IQ demodulation

         signal_Rx1 = RRC_filter2(I_Rx1,Q_Rx1,downsampling,beta,span); % RRC filtration
    
        % 2nd stream RRC filtration + IQ demodulation

        [I_Rx2,Q_Rx2] = IQ_demodulator(RF_signal_Tx2, carrier_real_Tx2, carrier_imag_Tx2); % IQ demodulation

        signal_Rx2 = RRC_filter2(I_Rx2,Q_Rx2,downsampling,beta,span); % RRC filtration

        
        %% Streams OFDM demodulator DATA + SIGNAL

        % 1st stream OFDM demodulation
        [ofdm_demodulated_Rx1, ofdm_demodulated_Rx1_pilots] = ofdm_demodulator2(signal_Rx1,transmission_type,CH_bandwidth,BCU_bandwidth,GI,Tx1_flag,equalizer_flag); % OFDM demodulation

        % 2nd stream OFDM demodulation
        [ofdm_demodulated_Rx2, ofdm_demodulated_Rx2_pilots] = ofdm_demodulator2(signal_Rx2,transmission_type,CH_bandwidth,BCU_bandwidth,GI,Tx2_flag,equalizer_flag); % OFDM demodulation

        ofdm_demodulated_output = [ofdm_demodulated_Rx1 ofdm_demodulated_Rx2];

        ofdm_demodulated_output_pilots = [ofdm_demodulated_Rx1_pilots ofdm_demodulated_Rx2_pilots];
        
        %% STBC 2Txx2Rx Alamouti code
        %ofdm_demodulated_Rx1 = Tx1;
        %ofdm_demodulated_Rx2 = Tx2;
        Rx = zeros(length(ofdm_demodulated_Rx2),2); % Rx streams

        for i = 1:length(ofdm_demodulated_Rx1)

             Rx(i,1) = ofdm_demodulated_Rx1(i); % 1st Tx stream
             Rx(i,2) = ofdm_demodulated_Rx2(i); % 2nd Tx stream
     %       Rx(i,1) = Tx1(i); % 1st Tx stream
     %       Rx(i,2) = Tx2(i); % 2nd Tx stream        
        end

        modulated_Rx_length = size(Rx,1); % length of modulated Rx stream

        modulated_Rx = []; % modulated Rx stream

        for i=1:2:modulated_Rx_length

            s1=Rx(i,1);
            s2=Rx(i+1,1);

            %Recieved data by RX1 Antenna at time interval T
            r(1,1)= (h(1,1)*s1) + (h(1,2)*s2); % + e(1,1);

            %Recieved data by RX1 Antenna at time interval (T+1)
            r(1,2)= ((-h(1,1))*conj(s2)) + (h(1,2)*conj(s1)); % + e(1,2);

            %Recieved data by RX2 Antenna at time interval T
            r(2,1)= (h(2,1)*s1) + (h(2,2)*s2); % + e(2,1);

            %Recieved data by RX2 Antenna at time interval (T+1)
            r(2,2)= ((-h(2,1))*conj(s2)) + (h(2,2)*conj(s1)); % + e(2,2);

            t(1,1)=((conj(h(1,1))*r(1,1)));
            t(1,2)=h(1,2)*(conj(r(1,2)));
            t(2,1)=((conj(h(2,1)))*r(2,1));
            t(2,2)=((h(1,2)*(conj(r(2,2)))));

            % channel matrix aproximation
            c(1,1)= ((conj(h(1,2)))*r(1,1));
            c(1,2)= h(1,1)*(conj(r(1,2)));
            c(2,1)= ((conj(h(2,2)))*r(2,1));
            c(2,2)= ((h(2,1)*(conj(r(2,2)))));

            %Maximum Likelehhod Detection Scehme
            s1_e =t(1,1) + t(1,2) + t(2,1) + t(2,2);
            s2_e= c(1,1) - c(1,2) + c(2,1) - c(2,2);

            Rx_per_period = [s1_e, s2_e];

            modulated_Rx = [modulated_Rx Rx_per_period];

        end
        % power_coeff = 7.694; % power coefficient for correcting modulation
        % modulated_Rx = modulated_Rx * power_coeff; 

    %     N_CBPS_signal = 48; % number of coded bits per symbol in BPSK modulation
    %     
    %     modulated_Tx1 = [modulated_Tx1_SIGNAL modulated_Tx1_DATA]; % 1st modulated Tx stream
    %     modulated_Tx2 = modulated_Tx2_DATA; % 2nd modulated Tx stream
    %     modulated_Tx = [modulated_Tx1 modulated_Tx2]; % modulated Tx stream 
    %     %power_coeff = 7.694; % power coefficient for correcting modulation
    %     %modulated_Tx = modulated_Tx * power_coeff;
    %     modulated_Tx_length = length(modulated_Tx); % length of modulated Tx stream

        N_CBPS = 48; % number of coded bits per symbol in BPSK modulation
        modulated_Rx1_SIGNAL = modulated_Rx(1:N_CBPS); % 1st modulated Rx stream (SIGNAL)
        modulated_Rx_DATA = modulated_Rx(N_CBPS+1:end); % Rx stream DATA
        modulated_Rx_DATA_length = length(modulated_Rx_DATA); % length of Rx stream DATA

        N_BPSC = T(MCS_index,:).N_BPSC; % read number of coded bits per subcarrier from T
        N_CBPS = T(MCS_index,:).N_CBPS; % read number of coded bits per symbol from T
        coded_DATA_Rx_length = modulated_Rx_DATA_length * N_BPSC; % number of coded Rx DATA bits
        number_of_blocks = coded_DATA_Rx_length / N_CBPS; % number of blocks

        % Recognition of DATA blocks contribution between streams

        coded_Rx1 = modulated_Rx1_SIGNAL; % 1st coded Rx stream (SIGNAL + DATA) 
        coded_Rx2 = []; % 2nd coded Rx stream (DATA)

        coded_data_zeros = zeros(1,coded_DATA_Rx_length);

        for block_number = 1:number_of_blocks
            if ( mod( block_number , 2 ) == 1) % odd blocks goes to 2nd Tx stream
                coded_Rx2 = [coded_Rx2  coded_data_zeros( ( N_CBPS * ( block_number - 1 ) )+1 : ( N_CBPS * block_number  ) ) ];
            else % even blocks goes to 1st Tx stream
                coded_Rx1 = [coded_Rx1 coded_data_zeros( ( N_CBPS * ( block_number - 1 ) )+1 : ( N_CBPS * block_number  ) ) ];
            end
        end

        modulated_Rx1_SIGNAL_length = length(modulated_Rx1_SIGNAL); % length of SIGNAL
        modulated_Rx1_DATA_length = ( length(coded_Rx1) - modulated_Rx1_SIGNAL_length ) / N_BPSC;
        modulated_Rx2_DATA_length = length(coded_Rx2) / N_BPSC;

        modulated_Rx1_DATA = modulated_Rx_DATA( 1 : modulated_Rx1_DATA_length ); % 1st Rx stream (DATA)
        modulated_Rx2_DATA = modulated_Rx_DATA( modulated_Rx1_DATA_length + 1 : end ); % 2nd Rx stream (DATA)


        %% Streams subcarrier demodulator DATA + SIGNAL

        % 1st stream demodulation

        demodulated_Rx1_SIGNAL = demodulator(modulated_Rx1_SIGNAL,'BPSK'); % demodulated 1st Rx stream (SIGNAL)   
        demodulated_Rx1_DATA = demodulator(modulated_Rx1_DATA,Modulation); % demodulated 1st Rx stream (DATA)

        % 2nd stream demodulation

        demodulated_Rx2_DATA = demodulator(modulated_Rx2_DATA,Modulation); % demodulated 2nd Rx stream (DATA)


        %% Streams deinterleaving DATA + SIGNAL

        % deinterleaving 1st Rx stream

        deinterleaved_Rx1_SIGNAL = deinterleaver(demodulated_Rx1_SIGNAL, CH_bandwidth, transmission_type,'BPSK'); % deinterleaved 1st Rx stream (SIGNAL)
        deinterleaved_Rx1_DATA = deinterleaver(demodulated_Rx1_DATA, CH_bandwidth, transmission_type, Modulation); % deinterleaved 1st Rx stream (DATA) 

        % deinterleaving 2nd Rx stream

        deinterleaved_Rx2_DATA = deinterleaver(demodulated_Rx2_DATA, CH_bandwidth, transmission_type, Modulation); % deinterleaved 2nd Rx stream (DATA) 

        %deinterleaved_signal = deinterleaved_Rx1_SIGNAL;
        %deinterleaved_data = [deinterleaved_Rx1_DATA deinterleaved_Rx2_DATA];
        
        %% Steam parser

        N_CBPS = T(MCS_index,:).N_CBPS; % read number of coded bits per symbol from T
        coded_data_blocks_Rx1 = length(deinterleaved_Rx1_DATA) / N_CBPS; % number of coded data blocks Rx1
        coded_data_blocks_Rx2 = length(deinterleaved_Rx2_DATA) / N_CBPS; % number of coded data blocks Rx2
        coded_data_blocks = ( length(deinterleaved_Rx1_DATA) + length(deinterleaved_Rx2_DATA) ) / N_CBPS; % number of coded data blocks

        Rx1_DATA = reshape(deinterleaved_Rx1_DATA.',N_CBPS,coded_data_blocks_Rx1).';
        Rx2_DATA = reshape(deinterleaved_Rx2_DATA.',N_CBPS,coded_data_blocks_Rx2).';

        deinterleaved_data = []; % coded data
        j1 = 1; j2 = 1;

        for block_number = 1:coded_data_blocks
            if ( mod( block_number , 2 ) == 1) % odd blocks get from data stream Rx2
                deinterleaved_data = [deinterleaved_data Rx2_DATA(j1,:)];
                j1 = j1 + 1;
            else % even blocks get from data stream Rx1
                deinterleaved_data = [deinterleaved_data Rx1_DATA(j2,:)];
                j2 = j2 + 1;
            end
        end    

        deinterleaved_signal = deinterleaved_Rx1_SIGNAL;

        %% FEC decoding + puncturing DATA + SIGNAL

        decoded_signal = fec_decoder(deinterleaved_signal,'1/2'); % convolutionally decoded signal

        decoded_data = fec_decoder(deinterleaved_data,R); % convolutionally decoded data

        %% Descrambler DATA

        descramblered_data = descrambler(decoded_data,s,half_byte); % descrambler DATA
        
        %% BER calculation
        
%         % pre Viterbi decoder
%         pre_RxData = [deinterleaved_signal deinterleaved_data]; % received bits before Viterbi decoder
%         nErrors = biterr(pre_TxData,pre_RxData); % calculate the number of bit errors
%         numErrs = nErrors; % number of error bits
%         %numBits = length(pre_TxData); % number of transmitted bits
%         pre_BER_actual = numErrs / pre_numBits; % BER before Viterbi decoder calculated
%         pre_BER_vector = [pre_BER_vector pre_BER_actual]; % fill BER vector with calculated BER before Viterbi decoder on each SNR
%     
%         % post Viterbi decoder
%         post_RxData = [decoded_signal decoded_data]; % received bits after Viterbi decoder
%         nErrors = biterr(post_TxData,post_RxData); % calculate the number of bit errors
%         numErrs = nErrors; % number of error bits
%         %numBits = length(TxData); % number of transmitted bits
%         post_BER_actual = numErrs / post_numBits; % BER after Viterbi decoder calculated
%         post_BER_vector = [post_BER_vector post_BER_actual]; % fill BER vector with calculated BER after Viterbi decoder on each SNR
        

        TxData = [coded_signal coded_data]; % transmitted bits
        
        % pre Viterbi decoder
        RxData = [deinterleaved_signal deinterleaved_data]; % received bits before Viterbi decoder
        nErrors = biterr(pre_TxData,RxData); % calculate the number of bit errors
        numErrs = nErrors; % number of error bits
        %numBits = length(pre_TxData); % number of transmitted bits
        pre_BER_actual = numErrs / pre_numBits; % BER before Viterbi decoder calculated
        pre_BER_vector = [pre_BER_vector pre_BER_actual]; % fill BER vector with calculated BER before Viterbi decoder on each SNR
    
        % post Viterbi decoder
        RxData = [decoded_signal decoded_data]; % received bits after Viterbi decoder
        nErrors = biterr(post_TxData,RxData); % calculate the number of bit errors
        numErrs = nErrors; % number of error bits
        %numBits = length(TxData); % number of transmitted bits
        post_BER_actual = numErrs / post_numBits; % BER after Viterbi decoder calculated
        post_BER_vector = [post_BER_vector post_BER_actual]; % fill BER vector with calculated BER after Viterbi decoder on each SNR        
        
        perc = perc + perc_step;
        waitbar(perc/100,w,sprintf('%d%% Receiver processing...',perc));
        
    
        %% MER + EVM calculation

        MER_act = MER_calc(modulated_Tx, ofdm_demodulated_output); % MER calculation
        MER_vector = [MER_vector MER_act]; % MER vector accumulation

        EVM_act = EVM_calc(modulated_Tx, ofdm_demodulated_output); % EVM calculation
        EVM_vector = [EVM_vector EVM_act]; % EVM vector accumulation
        
    end
    modulated_output = modulated_Tx;
end

end


%% Generation 

    perc = 90;
    waitbar(perc/100,w,sprintf('%d%% Graphics calculation...',perc));
   
    %% Transmit spectrum    

%     PSD = pwelch(RF_signal,[],[],8192,(8e+6)*upsampling,'power');
%     %PSD = pwelch(RF_signal,[],[],8192,8);
%     if mod(length(PSD),2) == 1
%         PSD(end) = [];
%     end
%     
%     if ( handles.BW == 8 )
%         fc=868-12;
%         f_BW=(handles.BW)*(40/(handles.BW));
%     elseif ( handles.BW == 7 )
%         fc=883.5-12;
%         f_BW=(handles.BW)*(35/(handles.BW));
%     elseif ( handles.BW == 6 )
%         fc=883-12;
%         f_BW=(handles.BW)*(30/(handles.BW));
%     end
%            
%            
%     %f_BW=8*(40/8);%*5.625;%*6.15;
%     f = linspace(fc-f_BW/2,fc+f_BW/2,length(PSD));
%     
%     axes(handles.axes1)
%     
%     y = 10*log10(PSD/max(PSD));
%     
%     plot(f,y);
%     
%     fc = 868;
%     
%     if ( strcmp(selected_RF_model,'AWGN') )
%         %Min = 863; Max = 873;
%         Min = fc - 5; Max = fc + 5;
%     else % case of Rician, Rayleigh and PI channels
%         %Min=846; Max=864;
%         Min = fc - 22; Max= fc - 4;
%     end

    PSD = pwelch(ofdm_output,[],[],8192,8,'power');
    if mod(length(PSD),2) == 1
        PSD(end) = [];
    end
    
    fc=0;
    f_BW=handles.BW;%*(40/8);%*5.625;%*6.15;
    f = linspace(fc-f_BW/2,fc+f_BW/2,length(PSD));
    y = 10*log10(PSD/max(PSD));
    
    axes(handles.axes1)
    plot(f,y);
    
    Min=-(handles.BW/2); Max=(handles.BW/2);
    xlim([Min,Max]);
    set(gca,'XTick',[Min : 1 : Max]);
    xlabel('Frequency [MHz]')
    ylabel('Normalized Magnitude [dB]')
    %title('Baseband spectrum on the transmitter output');
    
    counts_above_20dB = length(find(y > -20));
    bandwidth_above_20dB=counts_above_20dB*((8*(40/8))/length(y));

    %% Constellation diagram
    
    axes(handles.axes4)
    plot(real(ofdm_demodulated_output),imag(ofdm_demodulated_output),'.','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','r')
    grid on
    ax = gca;
%     y_min=ax.XLim(1);
%     y_max=ax.XLim(2);
%     ax.YLim = [y_min y_max];
    %axis square
    xlabel('In-phase')
    ylabel('Quadrature')
    title('Constellation diagram')
    
    
    %% Received spectrum
    %if ( (antenna_mode == '1x1 SISO') || (antenna_mode == '2x2 MIMO') )
        
        %[ofdm_data_Rx, data_flag] = ofdm_modulator2(ofdm_demodulated_data,transmission_type,CH_bandwidth,BCU_bandwidth,GI); % DATA OFDM modulation
        %if ( strcmp(transmission_type,'non-HT') )
        %    [ofdm_signal_Rx, signal_flag] = ofdm_modulator2(ofdm_demodulated_signal,transmission_type,CH_bandwidth,BCU_bandwidth,GI); % SIGNAL OFDM modulation
        %    ofdm_modulated_Rx = [ofdm_signal_Rx ofdm_data_Rx]; % complete DATA + SIGNAL OFDM symbols
        %elseif ( strcmp(transmission_type,'VHT') )
        %    ofdm_modulated_Rx = ofdm_data_Rx; % complete DATA OFDM symbols
        %end
        %PSD2 = pwelch(ofdm_modulated_Rx,[],[],8192,8);
        PSD2 = pwelch(ofdm_demodulated_output_pilots,[],[],8192,8);
        if mod(length(PSD2),2) == 1
            PSD2(end) = [];
        end
        
        fc=0;
        %fc=868-12;
        f_BW=handles.BW;%*(40/8);%*5.625;%*6.15;
        f_Rx = linspace(fc-f_BW/2,fc+f_BW/2,length(PSD2));
        
        %f_Rx=[-4094:4095]*8/8192;
        y_Rx = 10*log10(PSD2/max(PSD2));
    
    
    %elseif (antenna_mode == '2x2 MIMO')
    %    ofdm_demodulated_output;
    %end
    
    %%
    
    
    perc = (100 * bandwidth_above_20dB) / 8;

end

% print out error rate results
if (result == 1) % SNR static
    set(handles.edit9,'String', sprintf( '%0.2e', pre_BER_vector ) );
    set(handles.edit10,'String', sprintf( '%0.2e', post_BER_vector ) );
    set(handles.edit11,'String', sprintf( '%0.2f', MER_vector ) );
    set(handles.edit12,'String', sprintf( '%0.2f', EVM_vector ) );
    set(handles.edit13,'String', num2str(SNR));
else % SNR dynamic
    set(handles.edit9,'String', sprintf( '%0.2e', pre_BER_vector(end) ) );
    set(handles.edit10,'String', sprintf( '%0.2e', post_BER_vector(end) ) );
    set(handles.edit11,'String', sprintf( '%0.2f', MER_vector(end) ) );
    set(handles.edit12,'String', sprintf( '%0.2f', EVM_vector(end) ) );
    set(handles.edit13,'String', num2str( SNR(end) ) );
end

handles = guidata(hObject);

handles.rx_f = f_Rx;
handles.rx_y = y_Rx;

handles.rx_mod = ofdm_demodulated_output;
handles.tx_mod = modulated_output;
handles.mer = MER_vector;
handles.evm = EVM_vector;

handles.pre_ber = pre_BER_vector; % pre BER vector
handles.post_ber = post_BER_vector; % post BER vector

handles.snr = SNR; % SNR vector

% To save results in the Excell file...
% % A=handles.pre_ber;
% % xlswrite('TestA.xls',A);
% % B=handles.post_ber;
% % xlswrite('TestB.xls',B);
% % C=handles.snr;

if (result == 0) % SNR dynamic
SNR_BERbefore_BERafter_MER = [handles.snr' handles.pre_ber' handles.post_ber' handles.mer']; %'
xlswrite('Results.xls',SNR_BERbefore_BERafter_MER);
end

handles.d = DATA; % DATA binary

handles.a = scrambled_data; % scrambled DATA 

%save('spectrum_PI_6MHz.mat','f_Rx','y_Rx')
%save('constellation_AWGN_64QAM.mat','ofdm_demodulated_output','modulated_output')

guidata(hObject, handles);
elapsedTime = toc % time estimation of simulation
perc = 100;
waitbar(perc/100,w,sprintf('%d%% Calculations completed...',perc));
n=2;
pause(n)
close(w)

set(handles.edit19,'String', num2str(round(elapsedTime*1000)));

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

updateParameters(hObject,eventdata, handles);

function updateParameters(hObject, eventdata, handles)

% read selected BCU bandwidth by user from popupmenu4
BCU_bandwidths = get(handles.popupmenu4,'string');
selected_BCU_bandwidth = get(handles.popupmenu4,'Value');
BCU_bandwidth = BCU_bandwidths{selected_BCU_bandwidth};

% read selected channel configuration by user from popupmenu5
CH_bandwidths = get(handles.popupmenu5,'string');
selected_CH_bandwidth = get(handles.popupmenu5,'Value');
CH_configuration = CH_bandwidths{selected_CH_bandwidth};

% Calculate channel bandwith according to selected BCU bandwidth and channel configuration
if ( BCU_bandwidth == '8 MHz') BCU = 1;
    handles.BCU_W = 144;
    handles.BW = 8;
elseif ( BCU_bandwidth == '7 MHz') BCU = 2;
    handles.BCU_W = 168;
    handles.BW = 7;
elseif ( BCU_bandwidth == '6 MHz') BCU = 3;
    handles.BCU_W = 144;
    handles.BW = 6;
end

switch CH_configuration;
case 'TVHT_W' % User selects TVHT_W
   CH_bandwidth = ['8 MHz';'7 MHz';'6 MHz'];
   handles.IFFT_size = 128;
   ChBW_coeff = 1;
case 'TVHT_2W' % User selects TVHT_2W
   CH_bandwidth = ['16 MHz';'14 MHz';'12 MHz'];
   handles.IFFT_size = 512;
   ChBW_coeff = 2;
case 'TVHT_W+W' % User selects TVHT_W+W
   CH_bandwidth = ['8+8 MHz';'7+7 MHz';'6+6 MHz'];
   handles.IFFT_size = 128;
   ChBW_coeff = 2;
case 'TVHT_4W' % User selects TVHT_4W
   CH_bandwidth = ['32 MHz';'28 MHz';'24 MHz'];
   handles.IFFT_size = 1024;
   ChBW_coeff = 4;
case 'TVHT_2W+2W' % User selects TVHT_2W+2W
   CH_bandwidth = ['16+16 MHz';'14+14 MHz';'12+12 MHz'];
   handles.IFFT_size = 512;
   ChBW_coeff = 4;
end

% read transmission format choice defined in uibuttongroup7
selected_radiobutton = get(handles.uibuttongroup7,'SelectedObject');
selected_format = get(selected_radiobutton,'string');
result = strcmp(selected_format,'Non-HT');
if (result == 1)
    transmission_type = 'non-HT';
    
    % read Modulation Coding Scheme parameters for duplicate non-HT mode
    T=readtable('MCS.txt','Format','%s %s %s %s %f %f %f %s');
    
    handles.N_SD = 96; % data subcarriers
    handles.N_SP = 8; % pilot subcarriers
    handles.N_ST = handles.N_SD + handles.N_SP; % 104 subcarriers in total per BCU
else
    transmission_type = 'VHT';
    
    % read Modulation Coding Scheme parameters for duplicate VHT mode
    T=readtable('MCS_TVHT.txt','Format','%s %s %s %f %f %f %s');
    
    handles.N_SD = 108; % data subcarriers
    handles.N_SP = 6; % pilot subcarriers
    handles.N_ST = handles.N_SD + handles.N_SP; % 114 subcarriers in total per BCU
end

T.Properties.RowNames = T.MCS'; % define MCS rows

% gets selected MCS index by user from popupmenu1
MCS_indexes = get(handles.popupmenu1,'string');
selected_MCS_index = get(handles.popupmenu1,'Value');
MCS_index = MCS_indexes{selected_MCS_index};

% read number of coded bits per symbol
N_BPSC = T(MCS_index,:).N_BPSC;
handles.M = 2^(N_BPSC);

% read modulation
Modulation = T(MCS_index,:).Modulation;
Modulation = Modulation{1}; % converts cell to string

% read coding rate
R = T(MCS_index,:).R;
R = R{1}; % converts cell to string
handles.crate = str2num(R);

% read number of coded bits per subcarrier
N_BPSC = T(MCS_index,:).N_BPSC;

% read number of coded bits per OFDM symbol
N_CBPS = T(MCS_index,:).N_CBPS;

% read number of data bits per OFDM symbol
N_DBPS = T(MCS_index,:).N_DBPS;

%% calculate data rate

delta_f = handles.BW / handles.BCU_W; % Subcarrier frequency spacing kHz
T_FFT = 1/delta_f; % OFDM symbol duration [micros]

% gets selected antenna system cofiguration by user from popupmenu3
configurations = get(handles.popupmenu3,'string');
selected_configuration = get(handles.popupmenu3,'Value');
antenna_mode = configurations{selected_configuration};

switch antenna_mode
    case '1x1 SISO' % User selects 1x1 SISO
        AM_coeff = 1;
    case '2x1 MISO' % User selects 2x1 MISO
        AM_coeff = 1;
    case '2x2 MIMO' % User selected 2x2 MIMO
        AM_coeff = 2;
end

% gets cyclic prefix choice defined in uibuttongroup4
selected_radiobutton = get(handles.uibuttongroup4,'SelectedObject');
selected_GI = get(selected_radiobutton,'string');
result = strcmp(selected_GI,'Normal');
if (result == 1)
    T_GI = T_FFT/4; % GI duration [micros]
else
    T_GI = T_FFT/8; % GI duration [micros]
end
    
T_OFDM = T_FFT + T_GI; % full OFDM symbol period [micros]

data_rate = N_DBPS / T_OFDM; % data rate calculation
data_rate = round(data_rate); % round data rate on 1 decimal place
data_rate = ChBW_coeff * AM_coeff * data_rate; % antenna mode used (N_SS) or channel BW config
data_rate = num2str(data_rate); % conversion data rate to char

%textLabel = sprintf('Selected MCS index  %s Mbit/s', data_rate);
textLabel = sprintf('Channel bandwidth: %s\nChannel configuration: %s\nTransmission type: %s\nModulation: %s\nCoding rate: R=%s\nNumber of coding bits per subcarrier: %d\nNumber of coding bits per OFDM symbol: %d\nNumber of data bits per OFDM symbol: %d\nData rate: %s Mbps',BCU_bandwidth,CH_bandwidth(BCU,:),transmission_type,Modulation,R,N_BPSC,N_CBPS,N_DBPS,data_rate);
set(handles.text6, 'String', textLabel);

% Save the handles structure.
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hFig=figure;
%ButtonH=uicontrol('Parent',hFig,'Style','pushbutton','String','Close','Units','normalized','Position',[0.9 0.0 0.1 0.05],'Visible','on');
A=handles.rx_mod;
B=handles.tx_mod;
plot(real(handles.rx_mod),imag(handles.rx_mod),'.','MarkerSize',5,'MarkerEdgeColor','blue','MarkerFaceColor','r')
%scatterplot(modulated_data)
grid on
ax = gca;
y_min=ax.XLim(1);
y_max=ax.XLim(2);
ax.YLim = [y_min y_max];
%axis square
xlabel('In-phase')
ylabel('Quadrature')
title('Constellation diagram')
guidata(hObject, handles);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure;
subplot(2,1,1)
semilogy(handles.snr,handles.pre_ber,'LineWidth',1)
%hold on

M = handles.M; % Number of states modulation
k = log2(M); % Number of bits per symbol
coderate = handles.crate; % Code rate
K=7; % encoder length
G0=133; % generator polynomial G0
G1=171; % generator polynomial G1
t=poly2trellis(K,[G1 G0]); % define trellis
spect = distspec(t,1);
EbN0 = handles.snr - 10*log10(k) % Ratio of bit energy to noise power spectral density
%berAWGN = berawgn(EbN0,'psk',M,'nondiff')
if (M > 4)
    berCoded = bercoding(EbN0,'conv','hard',coderate,spect,'qam',M);
    berAWGN = berawgn(EbN0,'qam',M);
else
    berCoded = bercoding(EbN0,'conv','hard',coderate,spect,'psk',M,'nodiff');
    berAWGN = berawgn(EbN0,'psk',M,'nondiff');
end
%semilogy(handles.snr,berAWGN)
%semilogy(handles.snr,berAWGN)

grid on
xlabel('CNR [dB]')
ylabel('BER [-]')
title('BER before Viterbi decoding')
% A=handles.pre_ber;
% xlswrite('TestA.xls',A);
% B=handles.post_ber;s
% xlswrite('TestB.xls',B);
subplot(2,1,2)
semilogy(handles.snr,handles.post_ber,'LineWidth',1)
%hold on
%semilogy(handles.snr,berCoded)
grid on
xlabel('CNR [dB]')
ylabel('BER [-]')
title('BER after Viterbi decoding')

guidata(hObject, handles);

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure;
% stem(handles.a)
% grid on
% xlabel('Cislo bytu')
% ylabel('Hodnota bytu [dec]')
% title('DATA pole po skramblovani')

subplot(2,1,1)
%semilogy(handles.snr,handles.mer,'LineWidth',1)
plot(handles.snr,handles.mer,'LineWidth',1)
grid on
xlabel('CNR [dB]')
ylabel('MER [dB]')
title('Modulation error ratio')
C=handles.mer;
D=handles.evm;
subplot(2,1,2)
plot(handles.snr,handles.evm,'LineWidth',1)
grid on
xlabel('CNR [dB]')
ylabel('EVM [%]')
title('Error vector magnitude')

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function uibuttongroup4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uibuttongroup5.
function uibuttongroup5_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup5 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% gets choice  in uibuttongroup4
selected_radiobutton = get(handles.uibuttongroup5,'SelectedObject');
selected_SNR = get(selected_radiobutton,'string');
result = strcmp(selected_SNR,'Static');
if (result == 1)
    set(handles.edit7, 'Enable', 'on');
    set(handles.edit4, 'Enable', 'off');
    set(handles.edit5, 'Enable', 'off');
    set(handles.edit6, 'Enable', 'off');
else
    set(handles.edit7, 'Enable', 'off');
    set(handles.edit4, 'Enable', 'on');
    set(handles.edit5, 'Enable', 'on');
    set(handles.edit6, 'Enable', 'on');
end

function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
updateParameters(hObject,eventdata, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4
updateParameters(hObject,eventdata, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5
updateParameters(hObject,eventdata, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes when selected object is changed in uibuttongroup6.
function uibuttongroup6_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup6 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% gets value in radiobutton8 (AWGN channel)
% if active -> disable equalization checkbox1, otherwise enable it
value_radiobutton = get(handles.radiobutton8,'Value');
if (value_radiobutton == 1)
    set(handles.checkbox1, 'Enable', 'off');
else
    set(handles.checkbox1, 'Enable', 'on');
end

% --- Executes when selected object is changed in uibuttongroup7.
function uibuttongroup7_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup7 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% gets value in radiobutton12 (transmission type non-HT )
% if active -> change MCS indexes corresponding non-HT, otherwise change
% MCS indexes corresponding TVHT
value_radiobutton = get(handles.radiobutton12,'Value');

% gets selected MCS index by user from popupmenu1
MCS_indexes = get(handles.popupmenu1,'string');
selected_MCS_index = get(handles.popupmenu1,'Value');
%MCS_index = MCS_indexes{selected_MCS_index};

if (value_radiobutton == 1)
    MCS_indexes = num2cell(['0';'1';'2';'3';'4';'5';'6';'7']);
    set(handles.popupmenu1, 'String',MCS_indexes);
    if ( selected_MCS_index > 8 ) % if MCS index exceed range for selected format select upper limit value 
        set(handles.popupmenu1, 'Value', 8);
    else
        set(handles.popupmenu1, 'Value', selected_MCS_index);
    end
else
    MCS_indexes = num2cell(['0';'1';'2';'3';'4';'5';'6';'7';'8';'9']);
    set(handles.popupmenu1, 'String',MCS_indexes);
    set(handles.popupmenu1, 'Value', selected_MCS_index);
end

updateParameters(hObject,eventdata, handles);

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure;

%ofdm_demodulated_output = handles.rx_mod;
f = handles.rx_f;
y = handles.rx_y;
%[ofdm_data_Rx, data_flag] = ofdm_modulator2(ofdm_demodulated_output,transmission_type,CH_bandwidth,BCU_bandwidth,GI);
% PSD = pwelch(RF_signal,[],[],8192,8);
% if mod(length(PSD),2) == 1
%   PSD(end) = [];
% end
% 
% t=[-(length(PSD)/2):(length(PSD)/2)-1]*8/length(PSD);     
%     
% y=10*log10((PSD));
% 
% fc=868-12;
% f_BW=8*(40/8);%*5.625;%*6.15;
% f = linspace(fc-f_BW/2,fc+f_BW/2,length(PSD));
% 
% y = 10*log10(PSD/max(PSD));
%plot(f,y);

%%Frequency specifications:
%dF = Fs/N;                      % hertz
%f = -Fs/2:dF:Fs/2-dF;           % hertz
%%Plot the spectrum:

%plot(f,abs(X)/N);

%plot(real(ofdm_demodulated_output),'LineWidth',1);
plot(f,y);
%Min = -10; Max = 10;
%xlim([Min,Max]);
%set(gca,'XTick',[Min : 2 : Max]);
xlabel('Frequency [MHz]')
ylabel('Normalized Magnitude [dB]')
%title('Baseband spectrum on the OFDM demodulator output');

guidata(hObject, handles);


% --- Executes when selected object is changed in uibuttongroup4.
function uibuttongroup4_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup4 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateParameters(hObject,eventdata, handles);

function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'Enable', 'off');



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'Enable', 'off');


function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'Enable', 'off');


function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'Enable', 'off');


function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'Enable', 'off');


function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Close_Button.
function Close_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Close_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear all;
close all;
clc;
