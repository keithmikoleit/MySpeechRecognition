function [ mfcc, logmelcep, deltas, features ] = my_mfcc( speech, Tw, Ts, FFTL, Fs )
% returns mfcc for input speech
%   parameter speech = 16 kHz recorded speech
%   return mfcc = mfc coefficients for input recording


% Frame the signal into short frames.
% 25 msec frames * 16 kHz sample frequency = 400 samples per frame
% frame step of 10 msec = 160 samples
% For each frame calculate the periodogram estimate of the power spectrum.
% the periodogram is the power spectral estimate, which is the 
% magnitude of the FFT squared using a hamming window
window = Fs * Tw;
nooverlap = Fs * Ts;
nfft =(nextpow2(length(speech)));
L = length(speech);

%preemphasis coefficient
alpha = 0.97;

% Preemphasis filtering
speech = filter( [1 -alpha], 1, speech ); 
%fvtool( [1 -alpha], 1 );

% number of frames with last frame discard instead of padding
% because we know it will be silence
M = floor((L-window)/nooverlap+1);

% create hamming window
w = hamming(window);

%reserver space for frames
frames = zeros(window, M);
for idx = 0:M - 1
   frames(1:window,idx+1) = speech(nooverlap*idx+1:nooverlap*idx+window,1);
   % apply hamming window
   frames(1:window,idx+1)=frames(1:window,idx+1).*w;
end

% calculate the magnitude of the spectrum
MAG = abs( fft(frames,FFTL,1));

% Apply the mel filterbank to the power spectra, sum the energy in each filter.

% first we have to calculate the mel frequencies
upperf = 4000;
lowerf = 300;
nofilters = 26;

melfb = melfilterbank(lowerf, upperf, nofilters, FFTL, Fs);

% number of cols in MAG are actually the number of frequency bins
[Prows,bins] = size(MAG);
% number of columns in our mel filter bank is actually the number
% of filers we have
[melrows,melcols] = size(melfb);
% reserve spaces for calculations
melspectrum = zeros(melcols, bins);

% now multiply each PSD by each mel filter, and some up the resulting
% coefficients for each multiply.  This will give us one number that
% represents the energy in the bin
for idx = 1:bins
    for mel_idx = 1:nofilters
        % keep the first 256 values based on the filter size
        melspectrum(mel_idx,idx) = sum(MAG(1:melrows,idx).*melfb(:,mel_idx));
    end
end

% Take the logarithm of all filterbank energies.
mellogspectrum = log(melspectrum);
% replace inf values from log of 0 with 0s so we can take the dct
mellogspectrum(~isfinite(mellogspectrum))=0;

% reserve spaces for final features
nofeatures = 13;
mfcc = zeros(nofeatures, bins);
for idx = 1:bins
   % Take the DCT of the log filterbank energies.
   temp =  dct(mellogspectrum(:,idx));
   % Keep DCT coefficients 1-13, discard the rest.
   mfcc(1:nofeatures,idx) = temp(0:nofeatures - 1);
end

logmelcep = mellogspectrum;

% calculate delta features with linear least squared first-order estimator
deltas = zeros(nofeatures, bins);
M = 2;
for idx = 3:bins - 2
    for j = 1:nofeatures
        num = 0;
        den = 0;
        for t = -M:M
            num = num + t*mfcc(j, idx + t);
            den = den + t^2;
        end
        deltas(j, idx) = num / den;
    end
end

features = cat(1, mfcc, deltas);
% cut off the first 2 and last 2 bins due to lack of deltas
features = features(:,3:bins-2);

end

