function [normalisedVar]=LDA_normalise(data)
% normalisedVar = reshape(zscore(data(:)),size(data,1),size(data,2));
avgVar=nanmean(data);
sdVar=nanstd(data);
normalisedVar=(data-avgVar)./sdVar;
end 