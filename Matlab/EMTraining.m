% this script creates an HMM from L utterances of a word
L = 10;
N = 6; % estimated number of states in this word

% mu=[1,0.2, 0.3, 0.4, 0.8, 0.7, 0.9, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
% D=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
% mu(1:100)=0.5;
% D(1:100)=1;
% testmfcc = mvnrnd(mu,D,26);

% import L sound recordings for each word
% get audio data from wav files
% create MFCC matrices from each recording group these by word
s = cell(1,L);
features = cell(1, L);

Tw = 0.025; % 25 msec windows
Ts = 0.010; % 10 msec shift
FFTL = 512; % 512 point FFT

utterance = input('wavefile name', 's');
% get the base recording and MFCCs
[s{1},Fs] = audioread(sprintf('%s.wav', utterance));
% get the mfcc's
[amfcc, logmelcep, deltas, features{1}] = my_mfcc(s{1}, Tw, Ts, FFTL, Fs);
% get the rest of the same word recordings up to L utterances
for l = 2:L
%     test = sprintf('%s%d.wav', utterance, l);
    [s{l},Fs] = audioread(sprintf('%s%d.wav', utterance, l));
    [amfcc, logmelcep, deltas, features{l}] = my_mfcc(s{l}, Tw, Ts, FFTL, Fs);
end

% create first guess lambda parameters (5 total) for each word
    
    % create initial state probability
    Pi = zeros(N,1);
    % assume we start in the first state
    %Pi(1,1) = 1;
    Pi(1:N) = 1/N;
    
    % create state transition matrix
    A = zeros(N);
    %Initialize A with 0.5 Transition Values
    for i = 1:N - 1
        A(i,i) = 0.5;
        A(i,i+1) = 0.5;
    end
    %Last Transition has a value of 1
    A(N,N) = 1;
        
    %A = [.5 .5 0; 0 .5 .5; 0 0 1];                  %Transition Matrix is N by N
    
    % create initial emmision probabilty, which will be a multivariate
    % gaussian, to get this we need the mean vector and the diagonal
    % of the covariance matrix for the MFCC vectors
    
    %B = [.33 .33 .33 .33; .33 .33 .33 .33; .33 .33 .33 .33];%Emission Matrix is N by T
 
%     % we need one mean vector for each state in each utterance
%     u = cell(1,N,L);
%     for l = 1:L
%         for i = 1:N
%             u{1,i,l} = zeros(1,X); % reserve space for mean vector
%         end
%     end
%     
%     % intial estimate of mean vector
%     for l = 1:L
%         for q = 1:N
%             for x = 1:X % for each feature
%                 sum = 0;
%                 for t = 1:T
%                    sum = sum + features{l}(x,t); 
%                 end
%                 u{1,q,l}(1,x) = sum/T;
%             end
%         end
%     end
%      
%      % initial estimate of the covariance matrix
%     % we only need the diagonal, which is the variance, of the rows, so
%     % dimension = 2 normalized by the number of obervations
%     C = cell(1,N,L);
%     for l = 1:L
%         for i = 1:N
%             C{1,i,l} = var(features{l},0,2);  
%         end
%     end
    
    % calculate global mean and variance
    % global features
    gfeatures = [];

    Tw = 0.025; % 25 msec windows
    Ts = 0.010; % 10 msec shift
    FFTL = 512; % 512 point FFT

    files = dir('*.wav');
    for file = files'
        % get speech data for current file
        [temps,Fs] = audioread(file.name);
        % create MFCCs
        [amfcc, logmelcep, deltas, tempfeatures] = my_mfcc(temps, Tw, Ts, FFTL, Fs);
        % concatentate to global MFCCs
        gfeatures = cat(2,gfeatures, tempfeatures);
    end

    % X : number of features
    [X,T] = size(gfeatures); % number of features and frames

    % calculate global mean
    gu = zeros(X,1);

    for x = 1:X % for each feature
        sum = 0;
        for t = 1:T
           sum = sum + gfeatures(x,t); 
        end
        gu(x,1) = sum/T;
    end

    % caculalate global covariance
    gC = var(gfeatures,0,2);  

    % we need one mean vector and variance for each state in each utterance
    u = cell(1,N);
    for l = 1:L
        for i = 1:N
            u{1,i} = gu; % set to global mean
        end
    end
    
    C = cell(1,N);
    for i = 1:N
        C{1,i} = gC;  
    end
    
    % initial estimate of multivariate gaussian where the covariance
    % is just the diagonal and is positive semi-definite
    B = cell(1,L);
    for l = 1:L
        [X,T] = size(features{l});
        B{1,l} = zeros(N, T);
    end
    
    for l = 1:L
        for j = 1:N
            [X,T] = size(features{l});
            for t = 1:T
                mul = 1;
                sum = 0;
                for x = 1:X
                   mul = mul*C{j}(x);
                   sum = sum + (((features{l}(x,t) - u{j}(x))^2)/C{j}(x));
                end
%                 first = (1/((2*pi)^(X/2)*sqrt(mul)));
%                 second = exp(-(1/2)*sum);
                B{1,l}(j,t) = (1/((2*pi)^(X/2)*sqrt(mul)))*exp(-(1/2)*sum);
            end
        end
    end
    
%     % normalize gaussian values because they seem to come out too small
%     for l = 1:L
%         for t = 1:T
%             scalefactor = 0;
%             for j = 1:N
%                 scalefactor = scalefactor + B{1,l}(j,t);
%             end
%              for j = 1:N
%                 B{1,l}(j,t) = B{1,l}(j,t)/scalefactor;
%             end
%         end
%     end
        
% create alpha and beta matrices (1 combo per utterance)
alpham = cell(1,L);
betam = cell(1,L);
alphac = cell(1,L);
alpha = zeros(1,L);
beta = zeros(1,L);
graphingalpha = cell(1,L);
for l = 1:L
    [alpha(l),alpham{1,l}, nalpham, alphac{1,l}] = alpharec(A, B{1,l}, Pi);
    graphingalpha{l}(end + 1) =  alpha(1,l);
    [beta(l), betam{1,l}] = betarec(A, B{1,l}, Pi, alphac{1,l});
end

tao = 0.01; % certainty factor
iterations = 0;
% run the EM algorithm

% check the log ratio of HMM result, which is the log of
% the ration of the probability of x(1:T) from the current lambda
% to the previous.  We can get the p(x(1:T)) from our alpha
done = false;
while ~done
    
    % update lambda parameters for each word with update equations    
        % update the initial state probability
        % p(Q1 = j|lambda, x(1:T)) = gamma(j,t)
        % can't do this with function because eventually we will sum
        % up all utterances
        [N, T] = size(alpham{1,1}); % the number of states and observations in the HMM
        
        % reset Pi
        Pi = zeros(1,N);
%         for j = 1:N
%             for l = 1:L
%                 den = 0;
%                 for q = 1:N 
%                     den = den + alpham{1,l}(q,t)*betam{1,l}(q,t);
%                 end
%                 Pi(j) = Pi(j) + (alpham{1,l}(j,t)*betam{1,l}(j,t))/den;
%             end
%             Pi(j) = Pi(j) / L;
%         end
        for j = 1:N
            for l = 1:L
                Pi(j) = Pi(j) + HMMgammaposteri(j,1,alpham{1,l},betam{1,l});
            end
            Pi(j) = Pi(j) / L;
        end
        
        % update the state transition probability
        for i = 1:N
            for j = 1:N
                num = 0;
                den = 0;
                for l = 1:L
                    [N, T] = size(alpham{1,l});
                    for t = 2:T-1 % loop through T frames in utterance l
                        num = num + HMMxiposteri(i,j,t,alpham{1,l},betam{1,l},A,B{1,l},alphac{1,l});
                    end
                    for t = 2:T
                        den = den + HMMgammaposteri(i,t,alpham{1,l},betam{1,l});
                    end
                end
                A(i,j) = num/den;
            end
        end
        
        % update the u and C vectors, need to do this for each state
        for q = 1:N
            for x = 1:X
               num = 0;
               den = 0;
               for l = 1:L
                   [X,T] = size(features{l});
                   for t = 1:T
                      num = num + features{1,l}(x, t)*HMMgammaposteri(q, t, alpham{1,l},betam{1,l});
                      den = den + HMMgammaposteri(q, t, alpham{1,l},betam{1,l});
                   end
               end
               u{q}(x) = num / den;
            end
        end
        
%         for q = 1:N
%             for x = 1:X
%                num = 0;
%                den = 0;
%                for l = 1:L
%                    [X,T] = size(features{l});
%                    for t = 1:T
%                        num = num + ((features{1,l}(x, t) - u{q}(x))*(features{1,l}(x, t) - u{q}(x)))*HMMgammaposteri(q, t, alpham{1,l},betam{1,l});
%                        den = den + HMMgammaposteri(q, t, alpham{1,l},betam{1,l});
%                    end
%                end
%                C{q}(x, 1) = num/den;
%             end
%         end      
        for q = 1:N
%             for x = 1:X
               num = 0;
               den = 0;
               for l = 1:L
                   [X,T] = size(features{l});
                   for t = 1:T
                       num = num + ((features{1,l}(:, t) - u{q}(x))*(features{1,l}(:, t) - u{q}(x)).')*HMMgammaposteri(q, t, alpham{1,l},betam{1,l});
                       den = den + HMMgammaposteri(q, t, alpham{1,l},betam{1,l});
                   end
               end
               C{q}(:, 1) = diag(num/den);
%             end
        end  
        
        % generate a new emission probability matrix
        for l = 1:L
            for j = 1:N
                [X,T] = size(features{l});
                for t = 1:T
                    mul = 1;
                    sum = 0;
                    for x = 1:X
                       mul = mul*C{j}(x);
                       sum = sum + (((features{l}(x,t) - u{j}(x))^2)/C{j}(x));
                    end
    %                 first = (1/((2*pi)^(X/2)*sqrt(mul)));
    %                 second = exp(-(1/2)*sum);
                    B{1,l}(j,t) = (1/((2*pi)^(X/2)*sqrt(mul)))*exp(-(1/2)*sum);
                end
            end
        end
        
%         for j = 1:N
%             for t = 1:T
%                 mul = 1;
%                 sum = 0;
%                 for x = 1:X
%                    mul = mul*C{j}(x);
%                    sum = sum + (((features(x,t) - u{j}(x))^2)/C{j}(x));
%                 end
%                 first = (1/((2*pi)^(X/2)*sqrt(mul)));
%                 second = exp(-(1/2)*sum);
%                 B(j,t) = (1/((2*pi)^(X/2)*sqrt(mul)))*exp(-(1/2)*sum);
%             end
%         end
        
    % save the old probability value
    alphaold = alpha;

    % create new alpha and beta recursions from the lambda parameters
    for l = 1:L
        [alpha(l),alpham{1,l}, nalpham, alphac{1,l}] = alpharec(A, B{1,l}, Pi);
        graphingalpha{l}(end + 1) =  alpha(1,l);
        [beta(l), betam{1,l}] = betarec(A, B{1,l}, Pi, alphac{1,l});
    end
    
    loget = log(alpha/alphaold)
    iterations = iterations + 1
%     if abs(log(alpha/alphaold)) <= tao
%         done = true;
%     end
    if iterations > 30
        done = true
    end
%     else
%         % the previous values become the current values
%         alpha = alphanew;
%         alpham = alphamnew;
%         beta = betanew;
%         betam = betamnew; 
%     end

end

% we need to store and return the lambda parameters
% Pi - initial state matrix
% A - transition matrix
% u - mean vector per state
% C - variance vector per state 
% Note: B is calculated for u and C
save(sprintf('%slambda',utterance), 'A', 'Pi', 'u', 'C');

figure()
for l = 1:L
    hold on
    plot(graphingalpha{1,l})
end

disp()

sprintf('done with %s.wav', utterance)

    
