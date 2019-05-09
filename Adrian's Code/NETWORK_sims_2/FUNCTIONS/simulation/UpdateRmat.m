function out=UpdateRmat(sim,s)
%% function to update Rmat based on previous time iteration
% simple forward scheme of finite differences
% solving basic memristor equation
% s as in settings
% sim refers to data contained in previous simulation step (currents
% voltages, widths and resistances
out=struct();


Currents=sparse(triu(sim.Currents{1})'+triu(sim.Currents{1}));
OldR=sim.Rmat{1};
OldW=sim.Wmat{1};
PowFact=s.Pow;
[k,l,~]=find(triu(sim.AdjMat{1}));

%% update W
slope=s.Factor.*abs(Currents)-OldW./s.Tau;
parW=s.Step*slope;

newW=OldW+parW;

AW=newW-OldW;
alpha=sparse(k,l,s.Alpha*randn(length(k),1)+0.1*s.Alpha*ones(length(k),1),size(AW,1),size(AW,2));
alpha=alpha'+alpha;
sto_noise=alpha.*AW; 


newW=newW+sto_noise;

Pow=PowFact*(Currents.^2).*OldR;

newW=newW-Pow;


newW(newW>=s.MaxW)=s.MaxW;
newW(newW<=s.IniW & newW>0)=s.IniW;
newW(newW<0)=s.IniW;

newW=sparse(newW);

%% update R according to model
switch s.Model
    case 'HP' 
        % continuous linear resistances model
        fact=(newW-sim.AdjMat{1}.*s.IniW)./(s.MaxW-s.IniW);
        newR=sparse((s.Ron-s.Roff).*fact+s.Roff.*sim.AdjMat{1});       
    case 'Zdenka'
        %Discrete Resistances model
        newR=sim.Rmat{1};        
        newR(newW>=s.WForm)=s.Ron;
        newR(newR>=s.Ron & newW<=s.WDissolve)=s.Roff;
end

%% get conductance
[i,j,r]=find(newR);
newG=sparse(i,j,1./r);

%%out matrices
out.AdjMat=sim.AdjMat{1};
out.Gmat=newG;
out.Wmat=newW;
out.Rmat=newR;


end