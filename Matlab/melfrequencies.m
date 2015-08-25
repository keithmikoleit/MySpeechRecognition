% claculate mel frequencies for a given frequency range
% The formula for converting from frequency to Mel scale is:
% M(f) = 1125*ln(1 + f/700)
% To go from Mels back to frequency:
% iM (m) = 700*(e^(m/1125) - 1)
% to define the filter bank we find the upper and lower mels,
% cacluate the associated frequencies

% defaults are for a 16 kHz speech signal seperated into bins and
% psd's using a 512 point fft

upperf = 8000;
lowerf = 300;
upperm = 1125*log(1 + upperf/700);
lowerm = 1125*log(1 + lowerf/700);
nofilters = 10;
nfft = 512;
Fs = 16000;

m = linspace(lowerm, upperm, nofilters + 2);
h = 700*(exp(m/1125) - 1);

% round the frequencies to the nearest bins in the 512 point fft's
h = floor((nfft+1)*h/Fs)

% create an array of of filter banks
H = zeros(h(end)+1, nofilters);

% for each filter
for idx = 2:nofilters + 1
    % for each frequncy from 0 to the max frequency
    for f = 1:h(end)+1
        if f<h(idx-1)
            H(f, idx) = 0;
        elseif (f>=h(idx-1)&&f<=h(idx))
            H(f, idx) = (f-h(idx-1))/(h(idx)-h(idx-1));
        elseif (f>=h(idx)&&f<=h(idx+1))
            H(f, idx) = (h(idx+1)-f)/(h(idx+1)-h(idx));
        elseif f>h(idx+1)
            H(f, idx) = 0;    
        end
    end
end

plot(H)


