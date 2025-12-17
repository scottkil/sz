%% --- Read in raw data --- %
% clear all; close all; clc
filename = 'F:\Chrome Downloads\Sunrae_UCSF_ExampleData\ST0620250615_.csv'; % full file path
seizures = sz_findSeizures('filename',filename,'pband',[7 10],'ptCut',95); % auto detected of potential SWDs
ratioThresh = 1; % ratio for removing bad seizures 
seizures_better = sz_rmBadsz(filename, seizures, ratioThresh); % remove seizures not meeting this threshold