function [kmidx,k_centroid,mean_dist_within]=dist_weight_kmeans(inputs,k_seed,numb_k)



distmat=zeros(numb_k,length(inputs));
for i=1:numb_k
    distmat(i,:)=sqrt((inputs-inputs(k_seed(i))).^2);
    
end

[~,kmidx]=min(distmat,[],1); %mins along columns

for i=1:numb_k
    if ~isequal(sum(kmidx==i),0)
        k_centroid(i)=sum(inputs(kmidx==i))/sum(kmidx==i);
        mean_dist_within(i)=mean(sqrt((inputs(kmidx==i)-k_centroid(i)).^2));
    else
        k_centroid(i)=NaN;
        mean_dist_within(i)=NaN;
    end
end
mean_dist_within=mean(mean_dist_within);