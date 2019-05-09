function simout=KirchoffSolver(Sim,Electrodes,TimeInd)


%% EXTENDED KIRCHOFF METHOD
%% VOLTAGE SOurceS ARE TREATED AS KNOWN SOLUTIONS
%% BUT THE RANK OF THE MATRIX IS KEPT BY ADDING TWO CURRENTS(source,drain)
%% AS NEW VARIABLES IN THE EXTENDED MATRIX.

Numel=height(Electrodes);
OpenFlag=zeros(1,Numel);
for i=1:Numel
    OpenFlag(i)=Electrodes.OpenFlag{i}(TimeInd);
end
idx=find(OpenFlag==0);

%% fill extended matrix with electrodes
imat=Sim.Gmat;
imat=imat-diag(sum(imat,2));
ifix=sparse(length(imat),1);

for i=1:length(idx)
    ifix(end+1)=Electrodes.Value{idx(i)}(TimeInd);
    imat(end+1,Electrodes.PosIndex(idx(i)))=1;
    imat(Electrodes.PosIndex(idx(i)),end+1)=1;
end


%%solve

vsols=imat\ifix;

% retrieve currents

[rows,cols,~]=find(Sim.AdjMat);
Curr=Sim.Gmat;
voltdif=Sim.Gmat;
for i=1:length(rows)
    Curr(cols(i),rows(i))=(vsols(cols(i))-vsols(rows(i)))...
        *Curr(cols(i), rows(i));
    voltdif(cols(i),rows(i))=(vsols(cols(i))-vsols(rows(i)));
    
end

simout.SumRule=sum(sum(Curr,2));
k=1;
for i=1:Numel
    if isequal(OpenFlag(i),1)
        simout.(strcat('I',Electrodes.Name{i}))=NaN;
        simout.(strcat('V',Electrodes.Name{i}))=NaN;        
        continue;
    end
    simout.(strcat('I',Electrodes.Name{i}))=vsols(length(Curr)+k);
    simout.(strcat('V',Electrodes.Name{i}))=Electrodes.Value{i}(TimeInd);
    k=k+1;
end
simout.Currents=Curr;
simout.Voltages=vsols(1:end-length(idx));
simout.VoltDif=voltdif;




end