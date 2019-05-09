% 20/04/19
% Loading Function for Zdenkas's Code v0.1

function network = Load_Zdenka_Code()
%% Load Network Data:
%choose network to load
cd('D:\alon_\Research\POSTGRAD\PhD\CODE\Zdenka''s Code\atomic-switch-network-1.3-beta\asn\connectivity\connectivity_data');
waitfor(msgbox('Select the Network saved data'));
[FileName,PathName] = uigetfile('*.mat','Select the Network saved data');
f=fullfile(PathName,FileName);
network=load(f);
end 

