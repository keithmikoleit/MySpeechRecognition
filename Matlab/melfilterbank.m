function [ melfilters ] = melfilterbank( lowerf, upperf, nofilters, FFTL, Fs )
%melfilterbank returns a 2D array of mel filter's
%   lowerf : lower frequency bound
%   upperf : upper frequency bound
%   nofilters : no of filters in the bank
%   FFTL : length of FFT's the filter bank will be used with
%   Fs : Sampling frequency of signal the bank will be used with 

% claculate mel frequencies for a given frequency range
% The formula for converting from frequency to Mel scale is:
% M(f) = 1125*ln(1 + f/700)
% To go from Mels back to frequency:
% iM (m) = 700*(e^(m/1125) - 1)
% to define the filter bank we find the upper and lower mels,
% cacluate the associated frequencies

upperm = 1125*log(1 + upperf/700);
lowerm = 1125*log(1 + lowerf/700);

m = linspace(lowerm, upperm, nofilters + 2);
h = 700*(exp(m/1125) - 1);

% round the frequencies to the nearest bins in the FFTL point fft's
h = floor((FFTL+1)*h/Fs);

% create an array of of filter banks
melfilters = zeros(h(end)+1, nofilters);

% for each filter
for idx = 2:nofilters + 1
    % for each frequncy from 0 to the max frequency
    for f = 1:h(end)+1
        if f<h(idx-1)
            melfilters(f, idx - 1) = 0;
        elseif (f>=h(idx-1)&&f<=h(idx))
            melfilters(f, idx - 1) = (f-h(idx-1))/(h(idx)-h(idx-1));
        elseif (f>=h(idx)&&f<=h(idx+1))
            melfilters(f, idx - 1) = (h(idx+1)-f)/(h(idx+1)-h(idx));
        elseif f>h(idx+1)
            melfilters(f, idx - 1) = 0;    
        end
    end
end

end

