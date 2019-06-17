function [P, L,normalisedOutput,normalisedInput]=LDA_Apply(W, rawOutputs,rawInputs)
[normalisedOutput]=LDA_normalise_drain(rawOutputs);
[normalisedInput]=LDA_normalise_source(rawInputs);

L = [ones(size(normalisedOutput,1),1) normalisedOutput] * W';
P = exp(L) ./ repmat(sum(exp(L),2),[1 2]);
end

%NORMALISE ACROSS VALUES - take mean of the initial Inputs, SD -
%subtracting mean, dividing by SD - for each Drain seperately.

%OR 
% Normalise everything to a unit vector

%-- Look into LDA multi-class; otherwise All vs One classifier.