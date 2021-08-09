function joint_hist = jh_cap(data1,data2,resolution,xMin,xMax,yMin,yMax)
%JH       Joint Histogram
%
%   JOINT_HIST = jh_cap(data1,data2,resolution,xMin,xMax,yMin,yMax);
%
%	Creates a joint histogram from inputs 'data1' and 'data2' (must be
%	equal)
%
%   Inputs:
%       data1,
%           time series data organized as TIME x NODES.
%       data2,
%           time series data organized as TIME x NODES.
%       resolution,
%           amount of ticks on both axes
%       xMin,xMax,yMin,yMax,
%           cap for histogram
%
%   Outputs:
%       Joint Histogram,
%           
%


    ybins = linspace(yMax,yMin,resolution);
    xbins = linspace(xMin,xMax,resolution);
    xNumBins = numel(xbins);
    yNumBins = numel(ybins);

    Xi = round( interp1(xbins, 1:xNumBins, data1(:), 'linear', 'extrap') );
    Yi = round( interp1(ybins, 1:yNumBins, data2(:), 'linear', 'extrap') );
    Xi = max( min(Xi,xNumBins), 1);
    Yi = max( min(Yi,yNumBins), 1);
    joint_hist = accumarray([Yi(:) Xi(:)], 1, [yNumBins xNumBins]);

end



