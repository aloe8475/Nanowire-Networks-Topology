function [fp,sp]=psd_local(current,time)
        %time norm.
        time=time-time(1);
        Ts=time(2)-time(1);
        Fs =1/Ts; % Sampling frequency
        L = length(current); % Length of signal  
      
        %FFT
        xdft=fft(current);
        
        %norm factor
        
            %Norm=sqrt(L*Fs);
        Norm=max(time)*Fs;
             
        
        %%PSD 
        st=floor(length(current)/2+1);
        Pxx=1/(Norm)*abs(xdft(1:st,:)).^2;
        
        %out
        fp=0:Fs/L:Fs/2; %freq
        Pxx(2:end-1)=2.*Pxx(2:end-1); %freqs. are counted twice.
        sp=Pxx';
end