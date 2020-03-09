function seq=SplitAndAveKMat(kdat,div)

s=size(kdat,1);
if s<div
    div=s;
end
if div>length(kdat)
    div=1;
end
i=1;

seq=[];
while i<length(kdat)
    if i+div>length(kdat)
        break;
    end
     ave=mean(kdat(i:i+div-1,:),1);
     seq=[seq fuzzy_bool(minMax(ave,0.01),'labview')];
     i=i+div;
end
     ave=mean(kdat(i:end,:),1);
     seq=[seq fuzzy_bool(minMax(ave,0.01),'labview')];





end