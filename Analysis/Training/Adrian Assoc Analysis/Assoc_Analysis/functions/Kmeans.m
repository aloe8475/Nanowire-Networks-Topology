function [sequence,col,centroids]=Kmeans(kdataList,~,Kchoose,Knumb)
% 1 d kmeans
switch Kchoose
    case 'BestK'
        Kl=1:Knumb;
    case 'MaxK'
        Kl=Knumb;
end
sequence=[];div=10;
for i=1:length(kdataList)
    sequence=[sequence SplitAndAveKMat(kdataList{i},div)];
end

epochs=50;
kidxlist={};
centroidlist={};
mean_dist=[];
for i=1:length(Kl)
    for j=1:epochs
    k_seed=randi(length(sequence),1,Kl(i));
    [kmidx,k_centroid,mean_dist_within]=dist_weight_kmeans(sequence,k_seed,Kl(i));
    kidxlist{i,j}=kmidx;
    centroidlist{i,j}=k_centroid;
    meandist(i,j)=mean_dist_within;
    end
end


[~,I]=min(meandist(:));
[i,j]=ind2sub(size(meandist),I);
centroids=centroidlist{i,j};
kidx=kidxlist{i,j};

cols=lines(length(centroids));
col=zeros(length(sequence),3);
for i=1:length(sequence)
    col(i,:)=cols(kidx(i),:);
end




end