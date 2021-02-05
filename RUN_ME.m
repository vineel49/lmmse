% LMMSE based channel estimation for OFDM systems

close all
clear all
clc
SNR_dB = 40; % SNR per bit in dB
chan_len = 10; % number of channel taps
FFT_len = 256; % #subcarriers.
Pilot_num = 32; % number of pilots in each frame
CP_len = chan_len-1; % length of the cyclic prefix.
fade_var_1D = 0.5; % 1D fade variance of the channel impulse response
NUM_BIT = (FFT_len-Pilot_num)*2; % # data bits in each frame. i.e OFDM symbol.
NUM_FRAMES = 10^3; % simulations run (number of frames simulated)
%--------------------------------------------------------------------------
Pilot_Position = 1:FFT_len/Pilot_num:FFT_len; % pilot locations
Data_Position = 1:FFT_len; 
Data_Position(Pilot_Position) = []; % data location indices
%--------------------------------------------------------------------------
% SNR  PARAMETERS - OVERALL RATE IS 2
SNR = 10^(0.1*SNR_dB); % SNR in LINEAR SCALE
NOISE_VAR_1D = 0.5*2*2*chan_len*fade_var_1D/(2*SNR*FFT_len); % 1D AWGN NOISE VARIANCE 
NOISE_STD_DEV = sqrt(NOISE_VAR_1D); % NOISE STANDARD DEVIATION
%------------------- PRECOMPUTATIONS---------------------------------------
% cross correlation matrix.
Cross_Cor_Matrix = zeros(FFT_len-Pilot_num,Pilot_num); 
    for i1 = 1:FFT_len-Pilot_num
    for i2 = 1:Pilot_num
    Cross_Cor_Matrix(i1,i2) = Gen_autocorr(fade_var_1D,Data_Position(i1),Pilot_Position(i2),chan_len,FFT_len);
    end
    end
    
% autocorrelation matrix Rh_ls,h_ls
R_Input=zeros(Pilot_num,Pilot_num); % autocorrelation matrix.
    for i1=1:Pilot_num
    for i2=1:Pilot_num
       R_Input(i1,i2) = Gen_autocorr(fade_var_1D,Pilot_Position(i1),Pilot_Position(i2),chan_len,FFT_len);
    end
    end
 X = diag((1+1i)*ones(1,Pilot_num));   % all pilots are 1+1i.
    
 Auto_Cor_Matrix = R_Input + NOISE_VAR_1D*FFT_len*inv(X*X');    
%--------------------------------------------------------------------------
% Weight matrix    
Weight_Matrix = Cross_Cor_Matrix*inv(Auto_Cor_Matrix); 
%--------------------------------------------------------------------------
tic()
C_BER = 0; % bit errors in each frame
for FRAME_CNT = 1:NUM_FRAMES

% SOURCE
A = randi([0 1],1,NUM_BIT);

% QPSK MAPPING
QPSK_DATA_SIG = 1-2*A(1:2:end) + 1i*(1-2*A(2:2:end));

F_SIG_NO_CP = zeros(1,FFT_len); 
F_SIG_NO_CP(Data_Position) = QPSK_DATA_SIG;
F_SIG_NO_CP(Pilot_Position) = 1+1i; 

% IFFT 
T_SIG_NO_CP = ifft(F_SIG_NO_CP);

% INSERTING CYCLIC PREFIX
T_SIG_CP = [T_SIG_NO_CP(end-CP_len+1:end) T_SIG_NO_CP];

%---------------     CHANNEL      -----------------------------------------
% RAYLEIGH FREQUENCY SELECTIVE FADING CHANNEL
FADE_CHAN = normrnd(0,sqrt(fade_var_1D),1,chan_len)+1i*normrnd(0,sqrt(fade_var_1D),1,chan_len);

% AWGN
AWGN = normrnd(0,NOISE_STD_DEV,1,FFT_len+CP_len+chan_len-1)+1i*normrnd(0,NOISE_STD_DEV,1,FFT_len+CP_len+chan_len-1);

% channel output
T_REC_SIG = conv(T_SIG_CP,FADE_CHAN) + AWGN;

%----------------      RECEIVER  ------------------------------------------
% CP & TRANSIENT SAMPLES REMOVAL
T_REC_SIG(1:CP_len) = [];
T_REC_SIG_NO_CP = T_REC_SIG(1:FFT_len);

% PERFORMING THE FFT
F_REC_SIG_NO_CP = fft(T_REC_SIG_NO_CP);

% LMMSE BASED PILOT AIDED CHANNEL ESTIMATION (PACE)
Filter_Input = (F_REC_SIG_NO_CP(Pilot_Position))/(1+1i);
 
% Channel DFT at the data locations
F_H_Data_location = Weight_Matrix*(Filter_Input.');
F_H_Data_location = F_H_Data_location.'; % now a row vector

% ML DETECTION
QPSK_SYM = [1+1i 1-1i -1+1i -1-1i];
QPSK_SYM1 = QPSK_SYM(1)*ones(1,FFT_len-Pilot_num);
QPSK_SYM2 = QPSK_SYM(2)*ones(1,FFT_len-Pilot_num);
QPSK_SYM3 = QPSK_SYM(3)*ones(1,FFT_len-Pilot_num);
QPSK_SYM4 = QPSK_SYM(4)*ones(1,FFT_len-Pilot_num);

DIST = zeros(4,FFT_len-Pilot_num);
DIST(1,:)=abs(F_REC_SIG_NO_CP(Data_Position) - F_H_Data_location.*QPSK_SYM1).^2; 
DIST(2,:)=abs(F_REC_SIG_NO_CP(Data_Position) - F_H_Data_location.*QPSK_SYM2).^2;
DIST(3,:)=abs(F_REC_SIG_NO_CP(Data_Position) - F_H_Data_location.*QPSK_SYM3).^2;
DIST(4,:)=abs(F_REC_SIG_NO_CP(Data_Position) - F_H_Data_location.*QPSK_SYM4).^2; 

% COMPARING EUCLIDEAN DISTANCE
[~,INDICES] = min(DIST,[],1);
% MAPPING INDICES TO QPSK SYMBOLS
DEC_QPSK_MAP_SYM = QPSK_SYM(INDICES);
% DEMAPPING QPSK SYMBOLS TO BITS
DEC_A = zeros(1,NUM_BIT);
DEC_A(1:2:end) = real(DEC_QPSK_MAP_SYM)<0;
DEC_A(2:2:end) = imag(DEC_QPSK_MAP_SYM)<0;

% CALCULATING BIT ERRORS IN EACH FRAME
C_BER = C_BER + nnz(A-DEC_A);
end
toc()
% bit error rate
BER = C_BER/(NUM_BIT*NUM_FRAMES)
