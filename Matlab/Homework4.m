% record and audio signal
Fs = 16000;
y = audiorecorder(Fs,16,1); % creates the record object.
disp('start recording');
recordblocking(y,2); % records 4 seconds of speech
disp('stop recording');
s = getaudiodata(y);

sound(s, Fs)

% get audio data from wav file
%[s,Fs] = audioread(input('wavefile name', 's'));

figure()
plot(s)
xlabel( 'Sample' ); 
ylabel( 'Amplitude' );
title( 'Recorded Speech Waveform' );

C = 13;            % number of cepstral coefficients
D = 13;            % deltas

Tw = 0.025; % 25 msec windows
Ts = 0.010; % 10 msec shift
FFTL = 512; % 512 point FFT

% save the audio data as a wave file
audiowrite(input('wavefile name', 's'),s,Fs);

% get the mfcc's
[amfcc, logmelcep, deltas, features] = my_mfcc(s, Tw, Ts, FFTL, Fs);

% Plot cepstrum over time
figure('Position', [30 100 800 200], 'PaperPositionMode', 'auto', ... 
     'color', 'w', 'PaperOrientation', 'landscape', 'Visible', 'on' ); 

imagesc( [1:size(amfcc,2)], [0:(C+D)-1], amfcc ); 
axis( 'xy' );
xlabel( 'Frame index' ); 
ylabel( 'Cepstrum (0 - 12) and Delta (13) Index' );
title( 'Mel frequency cepstrum' );