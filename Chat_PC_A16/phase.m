function MF=phase(MF_signal, preamble_max_index, fsfd, preamble)
reduced_preamble1=MF_signal(preamble_max_index-10*fsfd:fsfd:preamble_max_index-4*fsfd);
realpart=find(real(reduced_preamble1) < 0);
reduced_preamble2=reduced_preamble1;
reduced_preamble2(realpart)=reduced_preamble1(realpart).*exp(1i*pi);
calc_angle=mean(angle(reduced_preamble2));
if(reduced_preamble1(1) < 0) 
	calc_angle = calc_angle+pi;    %jag får inte detta att funka, varför?
end
MF=MF_signal.*exp(-1i*calc_angle);

end