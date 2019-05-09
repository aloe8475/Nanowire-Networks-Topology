function [freqVector, freqSignal] = fourier_analysis(timeVector, signal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns the power spectrum density of a real signal. 
%
% ARGUMENTS: 
% time vector - should have an even length
% signal - to be transformed
%
% OUTPUT:
% freqVector
% freqSignal
%
% REQUIRES:
% none
%
% USAGE:
%{
    see tester_fourier_analysis
%}
%
% Authors:
% Ido Marcus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    L  = length(timeVector);                                % length of signal   (number)
    Fs = 1/(timeVector(2)-timeVector(1));                   % Sampling frequency (Hz)

    freqSignal = fftshift(fft(signal));                     % perform fft and shift
    freqSignal = abs(freqSignal/L);                         % Normalize and calculate power spectrum (no complex numbers)

    freqSignal = freqSignal(L/2:end);                       % given that the signal is real, one side is enough
    freqSignal(2:end) = 2*freqSignal(2:end);                % double it (only the zero frequency doesn't appear twice)

    freqVector = Fs*(0:(L/2))'/L;                            % construction of frequency axis  