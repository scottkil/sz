function seizures = sz_findSeizures(varargin)%pband, ptCut, ttv, eegChannel, targetFS)
%% sz_findSeizures Finds seizures in an EEG/LFP traces based on power thresholding
%
% INPUTS:
%   %----------------------Name-Value-Pairs----------------------%
%   Name: 'filename'     Value: full file path to file with EEG data
%   Name: 'pband'        Value: passband filter limits for seizure detection (default = [4 8])
%   Name: 'ttv'          Value: trough threshold value. This value is multipled by the standard
%                               deviation of the EEG to set a threshold that seizure troughs must
%                               pass (default = 3)
%   Name: 'ptCut'        Value: percentile cuttoff threshold. This value determines the
%                               threshold for detecting potential seizures by thresholding a power ratio (default = 95)
%   Name: 'eegChannel'   Value: channel number of the EEG (usually 1 or 2) (default = 1)
%   Name: 'tb'           Value: time buffer (in seconds) before and after seizures to grab  (default = 1)
%   %------------------------------------------------------------%
%
% OUTPUTS:
%   seizures - structure containing information about detected seizures
%
% Written by Scott Kilianski
% Updated 12/09/25
% ------------------------------------------------------------ %

%% Quick check to get the correct 'findpeaks' function because the chronux
% toolbox function of the same name sometimes shadows the MATLAB signal
% processing toolobox version (i.e. the version we actually want)
%-------------------------------------------------------------------------%
fpList = which('findpeaks.m','-all'); % find all 'findpeaks.m'
fpInd = find(contains(fpList,[filesep 'toolbox' filesep 'signal' filesep 'signal' filesep 'findpeaks.m']),1,'first'); %get location of the correct 'findpeaks.m'
currdir = cd; % save current directory path
cd(fileparts(fpList{fpInd})); % change directory to locatin of correct 'findpeaks.m'
fpFun = @findpeaks; % set function handle to correct 'findpeaks.m'
cd(currdir); % cd back to original directory

%-------------------------------------------------------------------------%
%% Parse inputs
validScalarNum = @(x) isnumeric(x) && isscalar(x);
default_filename = [];
default_pband = [4 8];
default_ptCut = 95;
default_ttv = 5;
default_eegChannel = 1;
default_targetFS = 200;
default_plotFlag = 1;
default_tb = 0.5; % time buffer (in seconds) - time to grab before and after each detected seizure
p = inputParser;
addParameter(p,'filename',default_filename,@(x) isstring(x) || ischar(x));
addParameter(p,'pband',default_pband,@(x) numel(x)==2);
addParameter(p,'ptCut',default_ptCut,validScalarNum);
addParameter(p,'ttv',default_ttv,validScalarNum);
addParameter(p,'eegChannel',default_eegChannel,validScalarNum);
addParameter(p,'targetFS',default_targetFS);
addParameter(p,'plotFlag',default_plotFlag);
addParameter(p,'tb',default_tb);
parse(p,varargin{:});
filename = p.Results.filename;
pband = p.Results.pband;
ptCut = p.Results.ptCut;
ttv = p.Results.ttv;
eegChannel = p.Results.eegChannel;
targetFS = p.Results.targetFS;
plotFlag = p.Results.plotFlag;
tb = p.Results.tb;
detectionParameters(1,:) = {'pband','ptCut','ttv','eegChannel'};
detectionParameters(2,:) = {pband,ptCut,ttv,eegChannel};

%% Load in data
fprintf('Loading data...\n')

if isempty(filename)
    [fn,fp,rv] = uigetfile({'*.csv'});
    if ~rv % if no file selected, end function early
        return
    else
        filename = fullfile(fp,fn);
    end
end

[fp, fn, fext] = fileparts(filename);           % get file name, path, and extension
if strcmp(fext,'.csv')
    x = readtable(fullfile(fp,[fn fext]));     % read in data from .csv
    EEG.data = table2array(x(:,eegChannel+1)); % store EEG/LFP data in proper format
    EEG.time = x.Time_s_;                       % store time data in proper format
    EEG.finalFS = round(1/diff(EEG.time(1:2))); % get sampling frequency (Hz)
else
    error('File type unrecognized. Use .csv file types only');
end
detectionParameters = [detectionParameters,{'finalFS';EEG.finalFS}];

%% Calculate spectrogram and threshold bandpower in band specificed by pband
fprintf('Calculating spectrogram and bandpower...\n');
frange = [0.5 50]; % Hz limits for spectogram
[f,t,specgram] = sz_spec(EEG,frange);% make spectrogram
R1 = f>=pband(1) & f<=pband(2);      % fundamental frequency range
pow1 = sum(specgram(R1,:),1);        % power over time in fundamental frequency range
R2 = f>=pband(1)*2 & f<=pband(2)*2;  % 2nd harmonic frequency range
pow2 = sum(specgram(R2,:),1);        % power over time in 2nd harmonic range
specRest = sum(specgram(~R1&~R2,1)); % power in the other frequencies of the spectrogram

sigPow = (pow1+pow2)./specRest;        % signal power (this is the power ratio used for detection)
sigPow_int = interp1(t,sigPow,...
    EEG.time,'linear','extrap');       % interpolate signal power to EEG timeframe
detThresh = prctile(sigPow_int,ptCut); % detection threshold
retThresh = detThresh/2;               % return-to-baseline threshold

% --- High pass filter at pband(1) and get instantaneous power --- %
[b, a] = butter(4, pband(1) / (EEG.finalFS/2), 'high'); % 4th order Butterworth highpass
eeg_filt = filtfilt(b, a, EEG.data); % apply filter

%% Find where power crosses threhold (determined by rising over detThresh and falling under retThresh)
pzit = 0.15; % gap length under which to merge (seconds)
mszt = .5; % minimum seizure time duration (seconds)
minTimeSamples = mszt*EEG.finalFS;

fprintf('Finding potential events...\n');
ts = sz_FindEvents(sigPow_int,EEG.time,detThresh,retThresh,pzit);
ts(:,1) = ts(:,1)-round(EEG.finalFS*tb); % adding time buffer to event starts
ts(:,2) = ts(:,2)+round(EEG.finalFS*tb); % adding time buffer to event ends
% --- check for events outside of recording boundaries and remove if necessary --- %
ts(ts(:,1)<=0,:) = [];
ts(ts(:,2)>numel(EEG.time),:) = [];
% -------------------------------------------------------------------------------- %
ts(diff(ts,1,2)<minTimeSamples,:) = []; % remove short events
startEnd_interp = EEG.time(ts);         % get start and end times of events

%% Find putative events, detect troughs, and store everything in structure(sz)
goodTRint = 0.05; % minimum interval between SWD troughs (seconds)
ttv = -std(eeg_filt)*ttv; % calculate trough threshold value (standard deviation * user-defined multiplier)
outfn = sprintf('%s%s%s_seizures.mat',fp,'\',fn);   % name of the output file

for ii = 1:size(ts,1)
    eegInd = ts(ii,1):ts(ii,2);
    seizures(ii).time = EEG.time(eegInd); % find EEG.time-referenced.
    seizures(ii).EEG = EEG.data(eegInd);
    if strcmp(fext,'.bin')
        seizures(ii).idx = EEG.idx(eegInd);
    end
    seizures(ii).type = 'Unclassified';
    [trgh, locs] = findpeaks(-seizures(ii).EEG); % find troughs (negative peaks)
    locs(-trgh>ttv) = []; % remove those troughs that don't cross the threshold (ttv)
    trgh(-trgh>ttv) = []; % remove those troughs that don't cross the threshold (ttv)

    %if troughs are too close together, only take the maximum negative one
    trInt = diff(locs) * (1/EEG.finalFS); % time interval between troughs
    tooSmall = find(trInt < goodTRint,1,'first'); %
    while tooSmall % run iteratively until all the too-close troughs are removed
        if trgh(tooSmall) > trgh(tooSmall+1)
            locs(tooSmall+1) = [];
            trgh(tooSmall+1)=[];    %remove the smaller of the too troughs
        else
            locs(tooSmall) = [];
            trgh(tooSmall)=[]; %remove the smaller of the too troughs
        end
        trInt = diff(locs) * (1/EEG.finalFS); % time interval between troughs
        tooSmall = find(trInt < goodTRint,1,'first'); %
    end
    seizures(ii).trTimeInds = locs; seizures(ii).trVals = -trgh; % store trough time (indices) and values in sz structure
    seizures(ii).filename = outfn;
    seizures(ii).parameters = detectionParameters;
end

tmp ={seizures(:).trTimeInds};
badLog = cellfun(@isempty,tmp);
seizures(badLog) = [];
startEnd_interp(badLog,:) = [];

%% Plotting trace, thresholds, and identified putative seizures
plotFlag = 1;
if plotFlag % plotting option
    figure; ax(1) = subplot(311);
    plot(EEG.time, EEG.data,'k','LineWidth',2); title('EEG');
    hold on
    % plot(get(gca,'xlim'),[ttv,ttv],'b','linewidth',1.5); hold off;
    ax(2) = subplot(312);
    % plot(t,bands.broadLow,'k','linewidth',2);
    % plot(EEG.time,sig1,'k','linewidth',2);
    plot(EEG.time,sigPow_int,'k','linewidth',2);
    % title(sprintf('Smoothed >%dHz envelope power',pband(1)));
    hold on
    plot(get(gca,'xlim'),[retThresh,retThresh],'r','linewidth',1.5);
    plot(get(gca,'xlim'),[detThresh,detThresh],'r','linewidth',1.5); hold off;
    ax(3) = subplot(313);
    imagesc(t,f,specgram);
    clim([0 prctile(specgram(:),99)])
    set(ax(3),'YDir','normal')
    title('Spectrogram'); xlabel('Time (sec)'); ylabel('Frequency (Hz)');
    linkaxes(ax,'x');
    axes(ax(1)); hold on;
    yl = get(gca,'YLim');
    xlim([EEG.time(1) EEG.time(end)])
    for ii = 1:size(startEnd_interp,1)
        patch([startEnd_interp(ii,:),fliplr(startEnd_interp(ii,:))],...
            [yl(1),yl(1),yl(2),yl(2)],'g',...
            'EdgeColor','none','FaceAlpha',0.25);
    end
end % plotting option end

[fp, fn, fext] = fileparts(outfn);
try % try statement here because sometimes saving fails due to insufficient permissions
    save(outfn,'seizures'); % save output into same folder as filename
    fprintf('%s%s saved in %s\n',fn,fext,fp)
catch
    fprintf('%s%s could not be saved in %s, likely because of insufficient permissions\n',fn,fext,fp)
end

end %main function end
