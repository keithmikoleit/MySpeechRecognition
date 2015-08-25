function [ gamma ] = HMMgammaposteri( q, t, alpham, betam )
% HMMgammaposteri This function computes the gamma posteri for a given 
%                 HMM defined by its alpha and beta matrices.  The 
%                 gamma posteri is defined as p(Q(t) = i|x(1:T)
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
%   q      : the current state
%   t      : the current time or frame
%   alpham : the alpha recursion matrix for the HMM
%   betam  : the beta recursion matrix for the HMM
%
%   The following paramters are returned by the function
%   gamma : the gamma posteri of the HMM at time t for state q
%
%
%   the gamma posteri is calculated as
%    gamma =
%    [alpham(q, t)*betam(q, t)]/prob(x(1:T))
%    [alpham(q, t)*betam(q, t)]/sum(alpham(q, t)*(betam(q,t), for q = 1:N)

[N, T] = size(alpham); % the number of states and observations in the HMM
gamma = 0;
for j = 1:N % the current state
    gamma = gamma + alpham(j,t)*betam(j,t);
end
gamma = (alpham(q, t)*betam(q,t)) / gamma;

end

