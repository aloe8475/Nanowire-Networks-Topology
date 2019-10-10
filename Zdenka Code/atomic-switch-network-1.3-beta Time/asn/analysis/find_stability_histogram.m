function [hist_bin_centers, hist_values] = find_stability_histogram(signal,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function finds the histogram (normalized as probability density 
% function) of the duration of time intervals between changes in the
% signal.
%
% ARGUMENTS: 
% signal - the signal to be analyzed.
% dt - time difference between conscutive measurements.

% OUTPUT:
% hist_bin_centers - x axis of the histogram (bin centers) (dimension like 
%                    dt)
% hist_values - the values of the histogram.
%
% REQUIRES:
% none
%
% Authors:
% Ido Marcus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    fig = figure('visible','off');
    
    bin_edges = (1:3:300)*dt;
    stability_histogram = histogram(diff([0; find(signal(1:end-1)~=signal(2:end)); length(signal)])*dt, bin_edges, 'Normalization','pdf');
    hist_values = stability_histogram.Values;
    hist_bin_centers  = mean([stability_histogram.BinEdges(1:end-1);stability_histogram.BinEdges(2:end)]);
    
    close(fig);

end