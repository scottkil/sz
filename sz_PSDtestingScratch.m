%% --- Read in raw data --- %
% clear all; close all; clc
% filename = 'F:\Chrome Downloads\Sunrae_UCSF_ExampleData\ST0720250613_.csv';
filename = 'F:\Chrome Downloads\Sunrae_UCSF_ExampleData\ST0620250615_.csv';
seizures = sz_findSeizures('filename',filename,'pband',[7 10],'ptCut',95);
seizures_better = sz_rmBadsz(filename, seizures);
figure;
histogram(PBtoLowRatio);
k = 0;
Fig1 = figure;

%%
k = k+1;
sax(2) = subplot(212);
plot(f,PSD(k,:),'k');
ylabel('Power');
xlabel('Frequeny (Hz)')
title(sprintf('PB to Low Ratio:  %.3f',PBtoLowRatio(k)));
sax(1) = subplot(211);
plot(seizures(k).time,seizures(k).EEG,'k');
ylabel('Voltage (A.U.)');
xlabel('Time (seconds)');
title(sprintf('Seizure %d',k));

set(Fig1().Children,'FontSize',20)