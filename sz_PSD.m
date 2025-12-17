function [PSD, f] = sz_PSD(eeg, SWDInds)
%% sz_PSD
%
% INPUTS:
%   - eeg - structure with following fields:
%       data - voltage data
%       time - same length as data but time values (in seconds)
%   - SWDInds - # SWD x 2 matrix with [start, end] indices of SWDs
%
% OUTPUTS:
%   PSD - # SWD x #freqs matrix. Normalized power spectrum for each SWD
%   f - frequency vector

%% === Function Body Here === %%
Fs = 1/diff(eeg.time(1:2));

% --- Set up parameters for PSD --- %
lfc = 1; % low frequency cutoff
hfc = 50; % high frequency cutoff
winSec   = 1;  % window length (seconds)
wlen     = max(32, round(winSec*Fs));     % window length in samples
nfft     = 4*2^nextpow2(max(wlen,256));   % FFT size in samples

window = hamming(wlen,'periodic');
desired_overlap = 0.95;
noverlap = round(desired_overlap*wlen);

tBuff = Fs*0; % time buffer for retrieving data before/after SWD

PSD = []; % intialize PSD matrix
for szi = 1:size(SWDInds,1)
        tmpidx = SWDInds(szi,:);
        winidx = (tmpidx(1)-tBuff):(tmpidx(2)+tBuff);
        if ~winidx(1)
            winidx(1) = 1;
        end
        x  = eeg.data(winidx);
        
        % t = (1:length(winidx))-1;
        % t = t/Fs;
        % t  = EEG.time(:);
        % Fs = 1/median(diff(t));                 % sampling rate (Hz)

        % x  = x - mean(x,'omitnan');             % remove DC
        % x  = fillmissing(x,'linear');           % handle any NaNs
        %

        % - Spectrogram code here -- %
        % winSec   = 2;                         % window length (seconds)
        % wlen     = max(32, round(winSec*Fs));   % window length (samples)
        % w        = hamming(wlen,'periodic');
        % desired_overlap = 0.98;
        % overlap  = round(desired_overlap*wlen);             % overlap
        % nfft     = 4*2^nextpow2(max(wlen,256));   % FFT size
        % %
        % [S,F,T,P] = spectrogram(x, w, overlap, nfft, Fs, 'psd');  % P: power/Hz
        % PdB = 10*log10(P + eps);                % dB scale
        % %
        % figure;
        % % imagesc(T, F, PdB); axis xy tight
        % sax(1) = subplot(211);
        % plot(t,x,'k');
        % sax(2) = subplot(212);
        % imagesc(T, F, P); axis xy tight
        % xlabel('Time (s)'); ylabel('Frequency (Hz)')
        % title('EEG Spectrogram');
        % ylim([0 50])
        % linkaxes(sax,'x');
        % % colormap turbo     % use 'parula' if preferred
        % -------------------------------------------------------%
    if numel(x) < wlen % if the event is too short
        continue
    end

        [pxx,f] = pwelch(x, window, noverlap, nfft, Fs, 'psd');
        cutLog = f >= lfc & f <=hfc;
        pxx = pxx(cutLog);
        f = f(cutLog);

        % --- Compute total power
        % df = f(2) - f(1);
        % lfc_bp1 = 4.5; % low frequency cutoff
        % hfc_bp1 = 7;
        % lfc_bp2 = 10;
        % hfc_bp2 = 12;
        %
        % cutLog1 = f >= lfc_bp1 & f <= hfc_bp1;
        % pxx_band1 = pxx(cutLog1);
        %     cutLog2 = f >= lfc_bp2 & f <= hfc_bp2;
        % pxx_band2 = pxx(cutLog2);
        % total_power = sum(pxx * df);              % ~ Parseval-consistent
        % band_power1  = sum(pxx_band1 * df);
        % band_power2 = sum(pxx_band2*df);
        % rel_bp1 = band_power1 / total_power; % fraction of total power in [lfc, hfc]
        % rel_bp2 = band_power2/ total_power;
        % metr(szi,1:3) = [total_power/numel(x), rel_bp1, rel_bp2];

        PSD(szi,:) = pxx./sum(pxx);
end
