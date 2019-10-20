const=[(1+1i) (1-1i) (-1-1i) (-1+1i)]/sqrt(2);

fs = 48000;                                           % sampling frequency [Hz]
Rb = 300;                                             % bit rate [bit/sec]
N = 432;                                               % number of bits to transmit
fc = 3000;                               % carrier frequency
span = 6;                                             % span is half of the interval length
M = length(const);                                    % Number of symbols in the constellation
bpsymb = log2(M);                                     % Number of bits per symbol
fsymb = Rb/bpsymb;                                    % Symbol rate [symb/s]
Tsymb = 1/fsymb;                                      % symbol time
fsfd = fs/fsymb;                                      % Number of samples per symbol (choose fs such that fsfd is an integer )easier [samples/symb]
Tsamp = 1/fs;                                         % Sampling time
Rs=fsymb


% PREAMBLE
preamble = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];

pack=[0 1 1 1 0 1 0 0 0 1 1 0 0 1 0 1 0 1 1 1 0 0 1 1 0 1 1 1 0 1 0 0];

m_idx=bi2de(buffer(pack, bpsymb)','left-msb')'+1;        %bits to symbol index
symbol = const(m_idx);                                   %Look up symbols using the indices

% PREAMBLE


%tx_signal=symb2sig(symbols_sent, const);
symb_upsample=upsample(symbol, fsfd);             %upsample signal
preamble_upsample=upsample(preamble, fsfd);             %upsample signal
pre_symbols = [preamble_upsample symb_upsample];       % preamble and symbols added

[pulse, ~] = rtrcpuls(0.6,1/fsymb,fs,span);         %creates pulse 
signal = conv(pulse,pre_symbols);                  %creates pulsetrain

time_vector = (0:length(signal) - 1)*Tsamp;
tx_signal = 2*signal.*exp(2*(-1i)*pi*fc*time_vector); 

tx_signal = tx_signal/max(abs(tx_signal));          %nomralize

%receiver

span = 6;                                           
[pulse, ~] = rtrcpuls(0.6,1/Rs,fs,span);        
    preamble_upsample = upsample(preamble,fsfd);

    % prev_index gives index of last checked recording, then updates it to
    % the current index- which is now index of last non-checked recording.
    % Saves the prev_index before overwritten if the preamble is found

    % current_data is the data of the recorder to be checked for preample
    % and data
    current_data = tx_signal;
    current_data = current_data(1:end);

    % time is a vector, the same length as the current_data, used in the
    % exponent that is multiplied with the current_data in order to
    % demodulate it
    time = 1/fs*(0:length(current_data)-1);
    exp_sig = sqrt(2)*exp(2*1i*pi*fc*time);

    % demodulation and normalizing
    exp_sig_data = current_data.*exp_sig'; 
    exp_sig_data = exp_sig_data/max(abs(exp_sig_data));

        preamble_data = conv(preamble_upsample,pulse);                   % convolve preamble with pulse 
        preamble_conv = conv(fliplr(preamble_data), real(exp_sig_data)); % convolve preamble with the data (the real part)    
        preamble_conv = abs(preamble_conv)./sqrt(2);                     % normalizing
        [preamble_peak, ~] = max(real(preamble_conv));                   % finding highest peak of the convoluted data
preamble_peak
