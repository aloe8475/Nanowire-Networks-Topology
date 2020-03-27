function [sclist,col,centroidlist]=scorer(avelist,targetlist,ScoreAlg,Kchoose,Knumb)
sclist=[];centroidlist=[];col=[];
switch ScoreAlg
    case 'CrossCorrelation'
        sclist=cross_correlator(avelist,targetlist);
    case 'EuclideanDistance'
        sclist=Euclid_distance(avelist,targetlist);
    case 'Kmeans'
        [sclist,col,centroidlist]=Kmeans(avelist,targetlist,Kchoose,Knumb);
end
end