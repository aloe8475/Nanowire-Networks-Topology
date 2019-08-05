CurrentVec = zeros(size(Connectivity.EdgeList,2),1);
for i = size(Connectivity.EdgeList,2)
    CurrentVec(i) = Output.Currents(Connectivity.EdgeList(1,i), Connectivity.EdgeList(2,i));    
end