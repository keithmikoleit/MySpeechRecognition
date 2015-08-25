function [ Pi ] = UpdatePi( alpham, betam )
%UpdatePi UpdatePi returns and updated initial probability vector for an
%         HMM defined by its alpha and beta matrices.
%   
%         The update function comes from the Baum-Welch algorithm for
%         expectation maximization.  Using langrange multipliers to
%         find a maximum using derivatives, we end up with the following
%         update equation:
%
%         p(Q(1) = j|lambda^g, x(1:T) = alpham(j,1)*betam(j,1)/
%                                       sum(alpham(j,1)*betam(j,1), over j)
%
%         this is the same as our gamma posteri evaluated at t = 1,
%         so we can use that function

[N, T] = size(betam); % the number of states and observations in the HMM
Pi = zeros(N);

for j = 1:N
    Pi(j) = HMMgammaposteri(j, 1, alpham, betam);
end
    
end

