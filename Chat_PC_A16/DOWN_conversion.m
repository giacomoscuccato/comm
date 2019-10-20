function [signal_at_baseband] = DOWN_conversion(signal_at_radiofrequency,frequency_carrier, T_sampling)
% this function downconverts a radiofrequency signal in a baseband one


time_vector = (0:length(signal_at_radiofrequency) - 1)*T_sampling;

signal_at_baseband = sqrt(2)*signal_at_radiofrequency.*cos(2*pi*frequency_carrier.*time_vector);



end

