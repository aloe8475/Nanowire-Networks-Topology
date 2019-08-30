%% EVOLUTION v0.1
%  28/08/2019
%  
%%  Overview (To Do):
% - Load Premade Networks x N (Joel)
% - Write Associative Learning Code (Alon + Joel) 
% - Write Cross-Correlation scoring Code (Alon) - done
%    -- Multiple Patterns - overall Cross-Correlation score 
%                         - learning rate
% - Write Noisy Recombination natural selection code (Ruomin)
%    -- Generate another set of premade networks with parameters from Noisy
%       Recombination (Ruomin)
% - Submit to our new overlords, the AgNWs
%% Generate Networks:
%
%
%
%
%
%%  Associative Learning Task (Alon + HELP ME)
% - Select target pattern 
% (Which Electrodes are on and which are off - combinations))
numElec=4;
elecArray=1:numElec;
targets=nchoosek(1:numElec,2);
% - Select training pattern (typically the same as target)
training=targets;

%create logicl array of which channels are on and which are off
for i = 1:length(numElec)
    boolTargets(i,:)=ismember(elecArray,targets(i,:));
end 
% - Write code to place electrodes (source x 4, drain x 4) in specific locations of the network
%   -- Same location for electrodes - how to do this? 
%   -- Start by spreading electrodes evenly based on spatial distribution
%      (x & y coordinates)


%%  - TRAINING: Send DC signal through chosen pattern electrodes, but remove non-chosen
%   drain electrodes UNTIL THRESHOLD CURRENT IS REACHED (Joel) - do this 10 times

% Send DC Signal to Training Electrodes.
% NOTE: only place the selected training electrodes in sources

%% - TESTING: Send DC signal through chosen pattern electrodes, but don't
%   remove non-chosen drain electrodes.

%% Compute Cross-Correlation scores + learning rate
% -Take average current across all times for each channel
% 
TimeIndex=round(sliderT.Value); %find the rounded time

avecurr=zeros(1,numElec);

avecurr(SelChan)=mean(currData(1:TimeIndex,SelChan),1);

%sclist=zeros(1,length(avecurr));
for i=1:length(avecurr)
    I_ave=avecurr{i}; %average current at the ith channel
    I_tar=boolTargets(i,:); %logical of which channels were used;
    
    un=(I_ave-mean(I_ave));
    dos=(I_tar-mean(I_tar));

    u=sum(un.*dos);
    d=sqrt(sum(un.^2)*sum(dos.^2));
    
    sclist(i)=u/d;
end 
%% Noisy Recombination Code (Ruomin)
%
%
%
%
%
