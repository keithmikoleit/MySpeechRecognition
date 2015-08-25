function [ alpha, alpham, nalpham, alphac ] = alpharec( A, B, Pi )
% alpharec This function returns the alpha matrix for an HMM
%
%   The following parameters are used for annotation
%   N    : Numer of states in HMM
%   T    : Number of observations (in the HW it is written Mw)
%   Q    : States in the HMM
%   x    : observations in the HMM
%   x(i) : the observation for frame/time i
%   j    : current state
%   i    : previous state
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
%   alpham : a NxT matrix representing the alpha value given (j,t)
%            which is the proposed state and the current time/frame
%   alpha  : a scalar that represents the summation of all values
%          : in alpham for the HMM

%   First we initialize the first column of alpha with:
%      alpham(1,j) = p(j)*p(x(1)|j)
%   The we recursively move forward across the matrix with:
%      alpham(t,j) = p(x(t)|j)*sum(alpha(i)*p(j|i), from i=1:N)   

t = 1; % the current frame/time
j = 1; % the current state
[N, T] = size(B); % the number of states and observations in the HMM
alpham = zeros(N, T);
nalpham = zeros(N, T);

% to avoid underflow the alpha matrix values are all normalized to 1
% these normalization factors are also passed to the beta matrices
alphac = zeros(1, T);

% code adapted from http://www.shokhirev.com/nikolai/abc/alg/hmm/HMM_cpp.html

% initialize the first column of the alpha vector;
for j = 1:N
   alpham(j, t) = Pi(j)*B(j,t);
   nalpham(j, t) = Pi(j)*B(j,t);
   alphac(1, t) = alphac(1, t) + alpham(j, t);
end
for j = 1:N
    alpham(j, t) = alpham(j,t)/alphac(1,t);
end
% recursively move forward across the matrix
for t = 1:T-1
    for j = 1:N
        sum1 = 0;
        sum2 = 0;
        for i = 1:N
            sum1 = sum1+alpham(i,t)*A(i,j);
            sum2 = sum2 + nalpham(i,t)*A(i,j);
        end
        alpham(j,t+1) = sum1*B(j,t+1);
        nalpham(j,t+1) = sum2*B(j,t+1);
        alphac(1, t+1) = alphac(1, t+1) + alpham(j, t+1);
    end
    % apply scalar values to each state at time t
    for j = 1:N
        alpham(j, t+1) = alpham(j,t+1)/alphac(1,t+1);
    end
end

% termination
% t = T;
% alpha = 0;
% for j = 1:N
%     alpha = alpha + alpham(j,t);
% end
% for the scaled vector we can now use the product of the scalars
% to find our total probability of obesrving the given events
alpha = sum(log(alphac));

end

