% find the most probable word from a speech utterance from zero to four

clear all
% load lambda parameters for each word in the lexicon
W = 5; % number of words in the lexicon
lambda = cell(1,4,W);
Pi_index = 1;
A_index = 2;
u_index = 3;
C_index = 4;

% load zero
load('zerolambda.mat')
lambda{1,Pi_index,1} = Pi;
lambda{1,A_index,1} = A;
lambda{1,u_index,1} = u;
lambda{1,C_index,1} = C;

% load one
load('onelambda.mat')
lambda{1,Pi_index,2} = Pi;
lambda{1,A_index,2} = A;
lambda{1,u_index,2} = u;
lambda{1,C_index,2} = C;

% load two
load('twolambda.mat')
lambda{1,Pi_index,3} = Pi;
lambda{1,A_index,3} = A;
lambda{1,u_index,3} = u;
lambda{1,C_index,3} = C;

% load three
load('threelambda.mat')
lambda{1,Pi_index,4} = Pi;
lambda{1,A_index,4} = A;
lambda{1,u_index,4} = u;
lambda{1,C_index,4} = C;

% load four
load('fourlambda.mat')
lambda{1,Pi_index,5} = Pi;
lambda{1,A_index,5} = A;
lambda{1,u_index,5} = u;
lambda{1,C_index,5} = C;

Fs = 16000;
y = audiorecorder(Fs,16,1); % creates the record object.

done = false;
while(done == false)

    cmd = input('q: quit, r: record speech, l: load speech \r\n', 's');
    
    if(cmd == 'q')
        done = true;
    else
        if(cmd == 'l')
            % load from a previously recorded speech
            % file
            [s,Fs] = audioread(input('wavefile name', 's'));
        elseif(cmd == 'r')
            % record the speech and then trim off the beginning and ending silence
            Fs = 16000;
            disp('start recording');
            recordblocking(y,2); % records 4 seconds of speech
            disp('stop recording');
            s = getaudiodata(y);

            % we wannt 100 msec before and after the utterance recording
            margin_seconds = 0.100;
            % margin in samples
            margin = Fs * margin_seconds;
            % set the threshold for non silence at an amplitude of 0.05
            threshold = 0.3;
            % find the first sample that has an amplitude greater than threshold.
            first = find(abs(s)>threshold,1, 'first');
            % check to see if the recording wasnt just silence
            if(first <= margin)
                first = 2 * margin;
            end
            % move back margin samples from the first valid sample
            first = first - margin;
            % find the last sample that has an amplitude greater than threshold.
            last = find(abs(s)>threshold,1, 'last');
            if(last > length(s) - margin)
                last = length(s) - 2 * margin;
            end
            % move forward margin samples from the last valid sample
            last = last + margin;

            % create a new truncated speech string
            s = s(first:last);
        end
        if(~isempty(s))
            % get the MFCCs for the uknown utterance
            Tw = 0.025; % 25 msec windows
            Ts = 0.010; % 10 msec shift
            FFTL = 512; % 512 point FFT

            [amfcc, logmelcep, deltas, features] = my_mfcc(s, Tw, Ts, FFTL, Fs);

            % array to store the final score against each word
            score = zeros(1,W);

            % for each word in the lexicon
            for w = 1:W
                % generate a B matrix with the MFCCs from the unknown utterance and
                % the mean and variance from the lexicon words lambda parameters

                Pi = lambda{1,Pi_index,w};
                A = lambda{1,A_index,w};
                u = lambda{1,u_index,w};
                C = lambda{1,C_index,w};

                % X : number of MFCCs
                % T : number of frames in the unknown utterance
                [X,T] = size(features);
                % N : number of states in the current word
                [N,N] = size(lambda{1,A_index,w});
                B = zeros(N, T);

                 for j = 1:N
                    for t = 1:T
                        mul = 1;
                        sum = 0;
                        for x = 1:X
                           mul = mul*C{j}(x);
                           sum = sum + (((features(x,t) - u{j}(x))^2)/C{j}(x));
                        end
                        B(j,t) = (1/((2*pi)^(X/2)*sqrt(mul)))*exp(-(1/2)*sum);
                    end
                end

                % generate an alpha matrix
                 [alpha,alpham,nalpham,alphac] = alpharec(A, B, Pi);

                % store the alpha value
                score(w) = alpha;

            end

            % return the max alpha value as the most probable word
            result = find(score == max(score));

            disp(score);

            fprintf('You said: %d \r\n', (result - 1));
        else
            disp('Nothing recorded')
        end
    end
end

disp(' Finished ')