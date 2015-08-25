% this script truncates all .wav files in a folder so that they only have
% 100 msec of silence at the beginning and the end

input('Do you wish to proceed?', 's')

files = dir('*.wav');
for file = files'
    %name = load(file.name);

    % get audio data from wav file
    %name = input('wavefile name', 's');
    [s,Fs] = audioread(file.name);

    figure()
    plot(s)

    % we wannt 100 msec before and after the utterance recording
    margin_seconds = 0.100;
    % margin in samples
    margin = Fs * margin_seconds;
    % set the threshold for non silence at an amplitude of 0.05
    threshold = 0.3;
    % find the first sample that has an amplitude greater than threshold.
    first = find(abs(s)>threshold,1, 'first');
    % move back margin samples from the first valid sample
    first = first - margin;
    % find the last sample that has an amplitude greater than threshold.
    last = find(abs(s)>threshold,1, 'last');
    % move forward margin samples from the last valid sample
    last = last + margin;

    % create a new truncated speech string
    news = s(first:last);

    figure()
    plot(news)

    % save the new audio data as a wave file
    audiowrite(file.name,news,Fs);

end