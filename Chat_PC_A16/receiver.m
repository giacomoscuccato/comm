% RECEIVER 
%
% This is the receiver structure that you will have to complete.
% The function: receiver(fc) is a setup function for the receiver. Here,
% the audiorecorder object is initialized (see help audiorecorder or
% MATLAB's homepage for more information about the object).
% 
% The callback function audioTimerFcn() is a callback function that is
% triggered on a specified time interval (here it is determined in the
% setup function, by the variable time_value)
% 
% Your task is to extend this code to work in the project!
%%

function [audio_recorder] = receiver(fc)

fs = 14000; %sampling frequency the issue we've found is that the headphones/soundcards are having problems playing the samples leaving us with a clipped signal which obviously doesn't work
audio_recorder = audiorecorder(fs,24,1);% Creates the recorder

%attach callback function
time_value = 1; % how often the function should be called in seconds
set(audio_recorder,'TimerPeriod',time_value,'TimerFcn',@audioTimerFcn); % attach a function that should be called every second, the function that is called is specified below.

%ADD USER DATA FOR CALLBACK FUNCTION (DO NOT CHANGE THE NAMES OF THESE VARIABLES!)
audio_recorder.UserData.receive_complete = 0; % this is a flag that the while loop in the GUI will check
audio_recorder.UserData.pack  = []; %allocate for data package
audio_recorder.UserData.pwr_spect = []; %allocate for PSD
audio_recorder.UserData.const = []; %allocate for constellation
audio_recorder.UserData.eyed  = []; %allocate for eye diagram
audio_recorder.UserData.index = 1;       % and index used with time_value to create frames of data      
audio_recorder.UserData.fc = fc;              % carrier frequency we get from the GUI
audio_recorder.UserData.preamble_found = 0;   % preamble flag, set to 1 when preamble is found


record(audio_recorder); %start recording
end


% CALLBACK FUNCTION
% This function will be called every [time_value] seconds, where time_value
% is specified above. Note that, as stated in the project MEMO, all the
% fields: pwr_spect, eyed, const and pack need to be assigned if you want
% to get outputs in the GUI.

% So, let's see an example of where we have a pulse train as in Computer
% exercise 2 and let the GUI plot it. Note that we will just create 432
% random bits and hence, the GUI will not be able to decode the message but
% only to display the figures.
% Parameters in the example:
% f_s = 22050 [samples / second]
% R_s = 350 [symbols / second]
% fsfd = f_s/R_s [samples / symbol]
% a = 432 bits
% M = 4 (using QPSK as in computer exercise)

function audioTimerFcn(recObj, ~, ~)
const = [-3-3i,-3-1i,-3+3i,-3+1i,-1-3i,-1-1i,-1+3i,-1+1i,3-3i,3-1i,3+3i,3+1i,1-3i,1-1i,1+3i,1+1i]/3;
%amp_inner=1/sqrt(2);
%amp_outer=1;
%const = [(1 + 1i)*amp_inner, 1i*amp_outer, (1i-1)*amp_inner, -1*amp_outer, (1-1i)*amp_inner, 1*amp_outer, (-1-1i)*amp_inner -1i*amp_outer];

fs = 14000;                                             % sampling frequency, same as above                                             % bit rate [bit/sec]
N = 432;                                                % number of bits to transmit
fc = recObj.UserData.fc;                                % carrier frequency from above init
span = 6;                                               % the span for our rrc
fsymb = 200;                                            % Symbol rate [symb/s]
fsfd = fs/fsymb;                                        % Number of samples per symbol (choose fs such that fsfd is an integer )easier [samples/symb]
preamble = [ 1,1,1,1,1,-1,-1,1,1,-1,1,-1,1 ];           % our preamble used for detecting the signal, we use barkers version
M=log2(length(const));

[pulse, ~] = rtrcpuls(0.6,1/fsymb,fs,span);             %creates the pulse used for creating the signal and demodulating the signal
upsample_preamble = upsample(preamble, fsfd);           %upsampling the preamble

prev_index = recObj.UserData.index;                     %stores the current index to use incase we find a preamble(to go back so we get the whole signal)
recObj.UserData.index=recObj.CurrentSample;             %stores the current index

signal_noise = getaudiodata(recObj);                    %stores audiodata 
signal_noise = signal_noise(prev_index:end);            %take out a frame of audiodata to be analysed

time = 1/fs*(0:length(signal_noise)-1);                 
exp_sig = sqrt(2)*exp(2*1i*pi*fc*time);                 %equation for Downsampling


exp_sig_data = signal_noise.*exp_sig';                  %downsampling
exp_sig_data = exp_sig_data/max(abs(exp_sig_data));     %normalization

if recObj.UserData.preamble_found == 0                  %go in this loop if preamble isn't found
    MF = conv(upsample_preamble,pulse);                 %create a signal that looks like the preamble in the signal
    corr_re = conv(fliplr(MF), real(exp_sig_data));     %try to locate preamble
    corr_re = abs(corr_re)./sqrt(2);
    [max_value, ~] = max(real(corr_re));                %find the max of the correlated preamble and signal,(peak = start of actual data and end of preamble)
    disp(max_value)
    if max_value > 12                                   %threshold to surpass to deem something a peak or max
    set(recObj, 'Timerperiod', 2)                       %change time to account for whole signal and not just preamble
    recObj.UserData.index = prev_index;                 %start from the precious index so the whole fram is included in calculations
    recObj.UserData.preamble_found=1;                   %flag for found preamble
    end
elseif recObj.UserData.preamble_found == 1              %enter this loop if the preamble is found
    MF=conv(upsample_preamble,pulse);                   %create a matched filter        
    corr_re = conv(fliplr(MF), exp_sig_data);           %redo some calculations but this time we want the index for where the max peak is
    corr_re = abs(corr_re)/sqrt(2);
    [~, delay] = max(real(corr_re));                    %index for the max peak
    
    signal=conv(fliplr(pulse),conj(exp_sig_data));      %same thing as for the preamble but this time the goal is to find general peaks in signal
    signal=signal/max(abs(signal));                     %normalise the signal
    signal=phase(signal,delay,fsfd,preamble);           %phase shift the signal 
    eye = signal(delay:delay+((N/M)*fsfd)-1-fsfd/M);    %take out only the signal(removes preamble) used to create eyediagram
    signal = signal(delay:fsfd:delay+((N/M)*fsfd)-1-fsfd/M);    %take out only the signal(removes reamble) and fsfd points apart    

    %minimum distace 
	euc_dist = abs(repmat(signal,1,length(const))-repmat(const, length(signal), 1)).^2;    
	[~, i] = min(euc_dist, [], 2);                 

	i = de2bi(i'-1,'left-msb');            %convert to bits
	indtmp = i';
	bits = indtmp(:)';
set(recObj, 'Timerperiod', 1)             %change time again to look for new preamble
recObj.UserData.preamble_found=0;

%------------------------------------------------------------------------------
% HOW TO SAVE DATA FOR THE GUI
%   NOTE THAT THE EXAMPLE HERE IS ONLY USED TO SHOW HOW TO OUTPUT DATA
%------------------------------------------------------------------------------

% % Step 1: save the estimated bits
recObj.UserData.pack = bits;
% 
% % Step 2: save the sampled symbols
recObj.UserData.const = signal;
% 
% % Step 3: provide the matched filter output for the eye diagram
recObj.UserData.eyed.r = eye;
recObj.UserData.eyed.fsfd = fsfd;

% Step 4: Compute the PSD and save it. 
%!!!! NOTE !!!! the PSD should be computed on the BASE BAND signal BEFORE matched filtering
[pxx, f] = pwelch(eye,1024,768,1024, fs); % note that pwr_spect.f will be normalized frequencies
f = fftshift(f); %shift to be centered around fs
f(1:length(f)/2) = f(1:length(f)/2) - fs; % center to be around zero
p = fftshift(10*log10(pxx/max(pxx))); % shift, normalize and convert PSD to dB
recObj.UserData.pwr_spect.f = f;
recObj.UserData.pwr_spect.p = p;

% % In order to make the GUI look at the data, we need to set the
% % receive_complete flag equal to 1:
recObj.UserData.receive_complete = 1; 
   end 
end
