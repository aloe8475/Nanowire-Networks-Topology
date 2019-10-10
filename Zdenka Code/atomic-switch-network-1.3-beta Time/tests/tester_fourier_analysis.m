close all;
clear; 

Fs = 1000;            % Sampling frequency (Hz)
T = 1/Fs;             % Interval between successice samples) (sec)
L = 1000;             % Length of signal
t = (0:L-1)*T;        % Time vector
X = 10*sin(2*pi*20*t)+5*sin(2*pi*50*t);

[freqVector, freqSignal] = fourier_analysis(t,X);

subplot(2,1,1)
plot(t,X);
title('signal vs. time');
xlabel('time (sec)');
ylabel('signal (a.u)');

subplot(2,1,2)
plot(freqVector,freqSignal)
title('signal PSD');
xlabel('Frequency (Hz)');
ylabel('PSD');