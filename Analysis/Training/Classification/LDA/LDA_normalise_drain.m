function [normalisedVar]=LDA_normalise_drain(data)
% normalisedVar = reshape(zscore(data(:)),size(data,1),size(data,2));
avgVar=nanmean(data);
sdVar=nanstd(data);
normalisedVar=(data-avgVar)./sdVar;
% 
% %thresholding:
% %Transform data:
% a=abs(data*10000);
% M = movmean(a,51);
% 
% normalisedVar=a>M;

end 