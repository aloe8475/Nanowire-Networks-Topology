function [normalisedVar]=LDA_normalise_source(data)
% normalisedVar = reshape(zscore(data(:)),size(data,1),size(data,2));
avgVar=nanmean(data);
sdVar=nanstd(data);
normalisedVar=(data-avgVar)./sdVar;

%thresholding:
%Transform data:
% a=(data(:,1)*10000);
% b=(data(:,2)*10000);
% Ma = movmin(a,50,'omitnan');
% Mb = movmin(a,50,'omitnan');
% normalisedVar(:,1)=a>Ma;
% normalisedVar(:,2)=b>Mb;
end 