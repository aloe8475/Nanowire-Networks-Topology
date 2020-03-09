function normArr=minMax(inArr,threshold)

m=max(inArr);

inArr(inArr<=threshold)=0;
inArr(inArr>threshold)=inArr(inArr>threshold)/m;

normArr=inArr;



end