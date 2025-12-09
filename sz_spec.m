function [F,T,P] = sz_spec(EEG,frange)

x = EEG.data;
Fs = 1/diff(EEG.time(1:2));

%- Spectrogram code here -- %
desired_overlap = 0.90;   % FFT window overlap
winSec   = 5;                         % window length (seconds)
wlen     = max(32, round(winSec*Fs));   % window length (samples)
w        = hamming(wlen,'periodic');
overlap  = round(desired_overlap*wlen);             % overlap
nfft     = 4*2^nextpow2(max(wlen,256));   % FFT size
%
[S,F,T,P] = spectrogram(x, w, overlap, nfft, Fs, 'psd');  % P: power/Hz
T = T+EEG.time(1);
PdB = 10*log10(P + eps);                % dB scale
Flog = F>frange(1) &  F<frange(2);
F = F(Flog);
P = P(Flog,:); 
% Optional Plotting %
% figure;
% % imagesc(T, F, PdB); axis xy tight
% sax(1) = subplot(211);
% plot(t,x,'k');
% sax(2) = subplot(212);
% drawnow;
% imagesc(T, F, P); axis xy tight
% xlabel('Time (s)'); ylabel('Frequency (Hz)')
% title('EEG Spectrogram');
% ylim([0 50])
% linkaxes(sax,'x');
% colormap turbo     % use 'parula' if preferred
% clim([0 max(P(:))*.8])
end % function end