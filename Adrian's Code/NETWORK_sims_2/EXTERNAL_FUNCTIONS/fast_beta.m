function [beta]=fast_beta(freq,amp,cutpoint)


freq=freq(freq>0 & freq<=cutpoint);
amp=amp(freq>0 & freq<=cutpoint);


 p=polyfit(log10(freq(2:end)),log10(amp(2:end)),1);
 beta=p(1);
end