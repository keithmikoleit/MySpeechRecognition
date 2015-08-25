function [ xi ] = HMMxiposteri( i, j, t, alpham, betam, A, B, alphac )
% HMMxiposteri This function computes the xi posteri for a given 
%                 HMM defined by its alpha and beta matrices.  The 
%                 xi posteri is defined as 
%   p(Q(t-1), Q(t), x(1:T)) = 
%    p(x(t)|q(t))*B(Q(t))*p(Q(t)|Q(t-1))*A(Q(t-1)
%
%   The following parameters are used for annotation
%   N    : Numer of states in HMM
%   T    : Number of observations (in the HW it is written Mw)
%   Q    : States in the HMM
%   x    : observations in the HMM
%   x(i) : the observation for frame/time i
%
%   The following parameters are passed to the function
%   j      : current state index
%   i      : previous state index
%   t      : the current time or frame
%   alpham : the alpha recursion matrix for the HMM
%   betam  : the beta recursion matrix for the HMM
%   A : NxN transition matrix describing the pobality that the next
%       state is Q(j) given that the previous state is Q(i), or a(ij).
%   B : NxT Observation probability matrix describing the probaility that 
%       we get observation x(t) given that the current state is Q(j)
%
%   The following paramters are returned by the function
%   xi : the xi posteri of the HMM at time t for state q
%
%   the xi posteri is calculated as
%    xi = B(j,t)*betam(j,t)*A(i,j)*alpham(i,t-1)

[N, T] = size(alpham); % the number of states and observations in the HMM
% den = 0;
% for q = 1:N % the current state
%     den = den + alpham(q,t)*betam(q,t);
% end
% 
% xiold = (B(j,t+1)*betam(j,t+1)*A(i,j)*alpham(i,t)/alphac(1,t+1));
% 
% xiold = xiold / den;
% now that we have scaled 
xi = HMMgammaposteri(i,t,alpham,betam)*A(i,j)*B(j,t+1)*betam(j,t+1);
xi = xi/(betam(i,t)*alphac(1,t+1));

if(~isfinite(xi))
    xi = 0;
end

end

