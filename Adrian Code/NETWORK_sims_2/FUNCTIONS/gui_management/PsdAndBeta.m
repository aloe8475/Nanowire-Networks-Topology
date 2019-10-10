function [F,PSD,beta_legend]=PsdAndBeta(T,II,cutpoint)

it=min(size(T));
for i=1:it
    [F(:,i),PSD(:,i)]=psd_local(II(:,i),T(:,i));
    beta_legend{i}=strcat('beta: ',num2str(fast_beta(F(:,i),PSD(:,i),cutpoint)));
end
end