% This is a minimal set of instructions to profile a program
% and visualize the results using Matlab's built in tools.

% Check current profiler status
ProfStatus = profile('status');

if strcmp(ProfStatus.ProfilerStatus, 'off'),
	profile on;
end

% Program to profile 
addpath(genpath('../examples/'))
NanowireNetworkAtomicSwitchDCBias

% Store data into struct
profData = profile('info');

% Read report 
profview(0, profData)