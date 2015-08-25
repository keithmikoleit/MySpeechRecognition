function [ beta, betam ] = betarec( A, B, Pi, alphac )
% betarec This function returns the beta matrix for an HMM
%
%   The following parameters are used for annotation
%   N    : Numer of states in HMM
%   T    : Number of observations (in the HW it is written Mw)
%   Q    : States in the HMM
%   x    : observations in the HMM
%   x(i) : the observation for frame/time i
%   j    : state index, i + 1
%   i    : state index, j - 1
%
%   The following parameters are passed to the function
%   A : NxN transition matrix describing the pobality that the next
%       state is Q(j) given that the previous state is Q(i), or a(ij).
%   B : NxT Observation probability matrix describing the probaility that 
%       we get observation x(t) given that the current state is Q(j)
%   Pi: 1xN initial probability matrix.  This matrix is non negative and
%       sums to 1.  It gives the probability that we start in Pi(n).
%
%   The following paramters are returned by the function
%   betam : a NxT matrix representing the alpha value given (j,t)
%            which is the proposed state and the current time/frame
%   beta  : a scalar that represents the summation of all values
%          : in alpham for the HMM

%   First we initialize the last column of beta with:
%      betam(k,T) = 1
%   The we recursively move backward across the matrix with:
%      betam(i,t) = sum(betam(j, t+1)*p(j|i)*p(j,t+1), for j from 1:N)

[N, T] = size(B); % the number of states and observations in the HMM
t = T; % the current frame/time
j = 1; % the current state
betam = zeros(N, T);

% code adapted from http://www.shokhirev.com/nikolai/abc/alg/hmm/HMM_cpp.html

% initialize the first column of the alpha vector;
for j = 1:N
   betam(j, t) = 1;
end

% recursively move backward across the matrix
for t = T-1:-1:1
    for i = 1:N
        sum = 0;
        for j = 1:N
            sum = sum+betam(j,t+1)*A(i,j)*B(j, t+1);
        end
        % in order to prevent underflow we scale the beta parameters
        % by the associated alpha matrix normalization factors
        % normalize each value in the matrix by its normalization factor
        betam(i,t) = sum / alphac(1,t+1);
    end
end

% termination
t = 1;
beta = 0;
for j = 1:N
    beta = beta + Pi(j)*B(j,t)*betam(j,t);
end

end