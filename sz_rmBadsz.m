function [seizures] = sz_rmBadsz(filename, seizures,ratioThresh)
%% sz_rmBadsz_
%
% INPUTS:
%   filename - path to original .csv file
%   seizures - structure output from sz_findSeizures
%   ratioThresh - threshold for cutting out bad seizures. This threshold is
%                 applied to the ratio of max passband power to max power lower than
%                 lower limit of passband
%
% OUTPUTS:
%   seizures - same as input except with bad entries removed


eegChannel = 1;
    x = readtable(filename);     % read in data from .csv
    eeg.data = table2array(x(:,eegChannel+1)); % store EEG/LFP data in proper format
    eeg.time = x.Time_s_;                       % store time data in proper format
    eeg.finalFS = round(1/diff(eeg.time(1:2))); % get sampling frequency (Hz)
%% --- Get SWDInds --- %
SWDInds = [];
for szii = 1:numel(seizures)
    timm = seizures(szii).time;
    [~,timmIDX] = ismember(timm,eeg.time);
    SWDInds(szii,:) = [timmIDX(1),timmIDX(end)];
end

[PSD, f] = sz_PSD(eeg, SWDInds);
pbF = seizures(1).parameters{2,1}; % passband limits

% --- Get max PSD value for frequencies below lower limit of passband --- %
lowFidx = f <  pbF(1);% find all frequencies lower than low cut of pass band
lowPSD = PSD(:,lowFidx); %get the corresponding PSD
maxLowP = max(lowPSD,[],2); % find the maximum value in that PSD

% --- Get max PSD value for frequencies within passband --- %
pbFidx = f >= pbF(1) & f <=  pbF(2);% find all frequencies between them
pbpPSD = PSD(:,pbFidx);     % PSD within the passband
maxPbP = max(pbpPSD,[],2); % maximum power in passband

PBtoLowRatio = maxPbP./maxLowP;
rmLog = PBtoLowRatio < ratioThresh;
seizures(rmLog) = [];


%%
