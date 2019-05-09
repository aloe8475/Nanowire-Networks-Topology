function [meanVal, slopeVal] = pinkNoiseAnalysis(timeAxis, signal, timeStamps, plotResults)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function extracts the exponent of the 1/f pattern of the PSD
% independently for different segments of a signal. It is meant to test
% whether there is a trend, in which different states of the network
% (ON\OFF) translate to different noise patterns.
%
% ARGUMENTS: 
% timeAxis - the time axis.
% signal - the signal.
% timestamps - Nx2 matrix specifying N segments of the timeAxis.
% plotResults - output flag.
%
% OUTPUT:
% meanVal - an NX1 vector of the mean value of the signal in each of the
%           segments.
% slopeVal - an NX1 vector of the 1/f exponent of the PSD in each of the
%            segments.
%
% REQUIRES:
% fourier analysis
%
% USAGE:
% see tester_pinkNoiseAnalysis
%
% Authors:
% Ido Marcus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    meanVal  = zeros(size(timeStamps,1),1);
    slopeVal = zeros(size(timeStamps,1),1);

    for i = 1 : size(timeStamps,1)
        % Extract current section:
        currIndices = timeAxis >= timeStamps(i,1) & timeAxis < timeStamps(i,2);
        if mod(sum(currIndices),2) == 1;
            currIndices(find(currIndices,1)) = 0;
        end            
        
        currTimeAxis = timeAxis(currIndices);
        currSignal   = signal  (currIndices);
        % Fourier trnasform current section:
        [freqAxis, freqSignal] = fourier_analysis(currTimeAxis-min(currTimeAxis), currSignal);
        
        % Linear fit of PSD in log-log scale:
        freqSignal = freqSignal(freqAxis~=0);
        freqAxis   = freqAxis(freqAxis~=0)';
        coefficients = polyfit(log10(freqAxis),log10(freqSignal),1);
        
        % Save results:
        meanVal(i)  =  mean(currSignal);
        slopeVal(i) = -coefficients(1);
        
        % Plot results:
        if plotResults
            fitPSD = 10^coefficients(2)*freqAxis.^(coefficients(1));
            
            subplot(size(timeStamps,1),2,(i-1)*2+1)
            plot(currTimeAxis,currSignal);
            title(sprintf('Signal segment #%d', i));
            xlabel('Time (sec)');
            ylabel('Signal (a.u.)');
            xlim(timeStamps(i,:));
            
            subplot(size(timeStamps,1),2,(i-1)*2+2)
            loglog(freqAxis,freqSignal,'b');
            hold on;
            loglog(freqAxis,fitPSD,'r');
            title(sprintf('Signal segment #%d - PSD', i));
            xlabel('Frequency (Hz)');
            ylabel('PSD (au^2)');
            set(gca,'YTick', 10.^(-20:5:20));
        end
    end
    
    % Plot results:
    if plotResults
        figure;
        subplot(2,1,1);
        plot(timeAxis,signal);
        legendEntries = cell(size(timeStamps,1)+1,1);
        legendEntries{1} = 'Signal';
        hold all;
        for i = 1 : size(timeStamps,1)
            currIndices = timeAxis >= timeStamps(i,1) & timeAxis < timeStamps(i,2);
            plot(timeAxis(currIndices),signal(currIndices),'LineWidth',2);
            legendEntries{i+1} = sprintf('Segment #%d', i);
        end
        title('Signal vs. time')
        xlabel('Time (sec)');
        ylabel('Signal (au)');
        legend(legendEntries,'Location','west');
        grid on;
        
        subplot(2,1,2);
        semilogx(meanVal,slopeVal,'*');
        title('Mean signal vs. 1/f noise slope');
        xlabel('Mean (au)');
        ylabel('1/f noise slope (au^2*sec)');
        grid on;
    end