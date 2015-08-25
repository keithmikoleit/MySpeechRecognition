% claculate the global mean and covariance of all .wav file 
% speech data in the current folder

input('Do you wish to proceed?', 's')
% loop through all .wav files in folder

% global features
gfeatures = [];

Tw = 0.025; % 25 msec windows
Ts = 0.010; % 10 msec shift
FFTL = 512; % 512 point FFT

files = dir('*.wav');
for file = files'
    % get speech data for current file
    [s,Fs] = audioread(file.name);
    % create MFCCs
    [amfcc, logmelcep, deltas, features] = my_mfcc(s, Tw, Ts, FFTL, Fs);
    % concatentate to global MFCCs
    gfeatures = cat(2,gfeatures, features);
end

[X,T] = size(gfeatures); % number of features and frames

% calculate global mean
u = zeros(X,1);

for x = 1:X % for each feature
    sum = 0;
    for t = 1:T
       sum = sum + gfeatures(x,t); 
    end
    u(x,1) = sum/T;
end

% caculalate global covariance
C = var(gfeatures,0,2);  


