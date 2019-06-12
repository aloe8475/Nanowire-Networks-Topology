
%Choose which simulations to load:
computer=getenv('computername');
switch computer
    case 'W4PT80T2'
        dataPath='C:\Users\aloe8475\Documents\GitHub\CODE\Adrian''s Code\NETWORK_sims_2\Saved Networks\Simulations Only';
    case ''
        dataPath='/suphys/aloe8475/Documents/CODE/Adrian''s Code/NETWORK_sims_2/Saved Networks/Simulations Only';
    case 'LAPTOP-S1BV3HR7'
        dataPath='D:\alon_\Research\POSTGRAD\PhD\CODE\Adrian''s Code\NETWORK_sims_2\Saved Networks\Simulations Only';
end
cd(dataPath);
waitfor(msgbox('Select the Simulation data'))
[FileName,PathName] = uigetfile('*.mat','Select the Simulation data');
f=fullfile(PathName,FileName);
%Load Data
load(f);
SelSims=combine_sim_data(SelSims);


%Save new data:

save([dataPath '/COMBINED_' FileName],'SelSims');

function SelSims=combine_sim_data(SelSims)

%% This function takes 3 simulations from main_sims.m and combines them into one.
%% Please note - only use this function if you have selected "RESUME SELECTED SIMULATION" in main_sims.m

numSims=length(SelSims);

temp.Domain=SelSims{1}.SelDomain;
temp.SelLayout=SelSims{1}.SelLayout;
temp.Settings=SelSims{1}.Settings;
temp.Settings=rmfield(temp.Settings,'Time'); %remove these values so we can reset them with all 3 simulations info
temp.Name=SelSims{1}.Name;
temp.DateCreated=SelSims{1}.DateCreated;

for i = 1:numSims
    if i == 1
        temp.Time=[SelSims{i}.Time];
        temp.Data=[SelSims{i}.Data];
        temp.SimInfo.MaxI=[SelSims{i}.SimInfo.MaxI];
        temp.SimInfo.MaxW=[SelSims{i}.SimInfo.MaxW];
        temp.SimInfo.MaxV=[SelSims{i}.SimInfo.MaxV];
        temp.SimInfo.MinV=[SelSims{i}.SimInfo.MinV];
        temp.SimInfo.MinW=[SelSims{i}.SimInfo.MinW];
        temp.SimInfo.MinI=[SelSims{i}.SimInfo.MinI];
        temp.Electrodes=[SelSims{i}.Electrodes];
        temp.Electrodes=removevars(temp.Electrodes,{'Value', 'OpenFlag'}); %remove these values so we can reset them with all 3 simulations info
        tempValue(i,:)=[SelSims{i}.Electrodes.Value{:}];
        tempFlag(i,:)=[SelSims{i}.Electrodes.OpenFlag{:}];
    else
        temp.Time=[temp.Time temp.Time(end) + SelSims{i}.Time];
        temp.Data=[temp.Data; SelSims{i}.Data];
        tempValue=[tempValue SelSims{i}.Electrodes.Value{:}];
        tempFlag=[tempFlag SelSims{i}.Electrodes.OpenFlag{:}];
        temp.SimInfo.MaxI=[temp.SimInfo.MaxI; SelSims{i}.SimInfo.MaxI];
        temp.SimInfo.MinI=[temp.SimInfo.MinI; SelSims{i}.SimInfo.MinI];
        temp.SimInfo.MaxW=[temp.SimInfo.MaxW; SelSims{i}.SimInfo.MaxW];
        temp.SimInfo.MinW=[temp.SimInfo.MinW; SelSims{i}.SimInfo.MinW];
        temp.SimInfo.MaxV=[temp.SimInfo.MaxV; SelSims{i}.SimInfo.MaxV];
        temp.SimInfo.MinV=[temp.SimInfo.MinV; SelSims{i}.SimInfo.MinV];
    end
end 
    clear SelSims
    %Settings
    SelSims.Settings=temp.Settings;
    SelSims.Settings.Time=length(temp.Time)/100;
    %All of this is pretty much the same except maybe DUTY
    
    %Time
    SelSims.Time=temp.Time;
    %Electrodes:
    SelSims.Electrodes=temp.Electrodes;
    A=reshape(tempValue,height(temp.Electrodes),length(tempValue)/height(temp.Electrodes))';
    SelSims.Electrodes.Value=num2cell(A,1)';
    B=reshape(tempFlag,height(temp.Electrodes),length(tempFlag)/height(temp.Electrodes))';
    SelSims.Electrodes.OpenFlag=num2cell(B,1)';
    %SelDomain
    SelSims.SelDomain=temp.Domain; %this is the same across all of the same type of simulations
    %SelLayout
    SelSims.SelLayout=temp.SelLayout; %this is the same across all of the same type of simulations
    %LastW
    %Data
    SelSims.Data=temp.Data;
    %Name
    SelSims.Name=temp.Name;
    %DateCreated
    SelSims.DateCreated=temp.DateCreated;
    %SimInfo
    SelSims.SimInfo.MaxI=max(temp.SimInfo.MaxI);
    SelSims.SimInfo.MinI=min(temp.SimInfo.MinI);
    SelSims.SimInfo.MaxW=max(temp.SimInfo.MaxW);
    SelSims.SimInfo.MinW=min(temp.SimInfo.MinW);
    SelSims.SimInfo.MaxV=max(temp.SimInfo.MaxV);
    SelSims.SimInfo.MinV=min(temp.SimInfo.MinV);
    SelSims={SelSims};

end