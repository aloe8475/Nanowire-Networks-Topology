function [meanVal,varVal,binSize] = taylorLaw(signal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function implements the "expanding bins" method in order to examine
% whether a signal satisfies "Taylors Law", which is related to 1/f noise.
%
% ARGUMENTS: 
% signal - the input signal.
%
% OUTPUT:
% meanVal, varVal - the mean and variance values obtained.
%
% REQUIRES:
% none
%
% Authors:
% Ido Marcus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    N = length(signal);
    minSamples = 10;
    minBins    = 10;
    numPoints  = 100;
    binSize    = floor(logspace(log10(minSamples),log10(N/minBins),numPoints));

    meanVal = zeros(size(binSize));
    varVal  = zeros(size(binSize));

    for i = 1 : length(binSize)
        % First segment signal (every column is a segment), then sum:
        summedSegmentedSignal = sum(reshape(signal(1:end-mod(N,binSize(i))), [binSize(i),floor(N/binSize(i))])); 

        % Find statistics:
        meanVal(i) = mean(summedSegmentedSignal);
        varVal (i) = var (summedSegmentedSignal);
    end

%     loglog(meanVal,varVal,'*');
%     title('Variance vs. mean of signal');
%     xlabel('Mean');
%     ylabel('Variance');
%     grid on;
end