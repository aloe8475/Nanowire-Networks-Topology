%% TDMS Import
%choose TDMS file to import
computer=getenv('computername');
switch computer
    case 'W4PT80T2'
        dataPath='C:\Users\aloe8475\Documents\PhD\GitHub\CODE\Analysis\Training\Assoc_Analysis\figures\Training Data';
    case ''
        dataPath='/suphys/aloe8475/Documents/CODE/Analysis/Training/Assoc_Analysis/figures/Training Data';
    case 'LAPTOP-JTHFN66G'
        dataPath='Z:\Documents\CODE\Analysis\Training\Adrian Assoc Analysis\Assoc_Analysis\figures\Training Data';
end
cd(dataPath)
waitfor(msgbox('Select the TDMS saved data'));
[FileName,PathName] = uigetfile('*.mat','Select the TDMS saved data');
f=fullfile(PathName,FileName);

[ConvertedData,ConvertVer,ChanNames]=convertTDMS(true,f);