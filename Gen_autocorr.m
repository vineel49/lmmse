% Generates 1D autocorrelation of channel DFT in OFDM systems
% (37) in the paper K Vasudevan, “Coherent Detection of Turbo Coded OFDM Signals Transmitted through Frequency Selective Rayleigh Fading Channels ”, 
% IEEE International Conference on Signal Processing Computing and Control, 26—28 Sept. 2013, Shimla.
% 'l' and 'm' are the channel DFT location indices
function [auto_cor] = Gen_autocorr(fade_var_1D,l,m,chan_len,FFT_len)
auto_cor = 0; % autocorrelation initialization
for auto_cnt = 0:chan_len-1
   auto_cor = auto_cor + (fade_var_1D *exp(-1i*2*pi*(l-m)*auto_cnt/FFT_len));
end
end