function [normalisedVar]=LDA_normalise(normaliseVariable)
avgVar=nanmean(normaliseVariable);
sdVar=nanstd(normaliseVariable);
normalisedVar=(normaliseVariable-avgVar)./sdVar;
end 