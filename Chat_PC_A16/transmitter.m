% COMPLETE THE TRANSMITTER!

% pack = message to be transmitted (consists of 432 bits from the GUI, always!)
% fc = carrier frequency

function transmitter(pack,fc)
%amp_inner=1/sqrt(2);
%amp_outer=1;
%const = [(1 + 1i)*amp_inner, 1i*amp_outer, (1i-1)*amp_inner, -1*amp_outer, (1-1i)*amp_inner, 1*amp_outer, (-1-1i)*amp_inner -1i*amp_outer];
const = [-3-3i,-3-1i,-3+3i,-3+1i,-1-3i,-1-1i,-1+3i,-1+1i,3-3i,3-1i,3+3i,3+1i,1-3i,1-1i,1+3i,1+1i]/3;


fs = 14000;                                             % sampling frequency, same as above                                             % bit rate [bit/sec]
N = 432;                                                % number of bits to transmit
span = 6;                                               % the span for our rrc
fsymb = 200;                                            % Symbol rate [symb/s]
fsfd = fs/fsymb;                                        % Number of samples per symbol (choose fs such that fsfd is an integer )easier [samples/symb]
M=log2(length(const));
preamble = [ 1,1,1,1,1,-1,-1,1,1,-1,1,-1,1 ];           % our preamble used for detecting the signal, we use barkers version

[pulse, ~] = rtrcpuls(0.6,1/fsymb,fs,span);             %creates the pulse used for creating the signal and demodulating the signal

m_idx=bi2de(buffer(pack, M)','left-msb')'+1;            %bits to symbol
symbol = const(m_idx);                                  %Look up symbols using the indices

symb_upsample=upsample(symbol, fsfd);                   %upsample signal
preamble_upsample=upsample(preamble, fsfd);             %upsample preamble
pre_symbols = [preamble_upsample symb_upsample];        %put preamble first followed by signal

signal = conv(pulse,pre_symbols);                       %crete a pulsetrian(signal)

time_vector = (0:length(signal) - 1)*1/fs;              %time vector used for upconversion
tx_signal = 2*signal.*exp(2*(-1i)*pi*fc*time_vector);   %upconversion

tx_signal = tx_signal/max(abs(tx_signal));              %normalize the signal

player = audioplayer(real(tx_signal), fs);              %Audioobject that plays signal at fs frequency
playblocking(player);                                   %Play the signal



