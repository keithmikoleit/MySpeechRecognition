% test script for HMM

Pi = [.33 .33 .33];
A = [.5 .5 0; 0 .5 .5; 0 0 1];                  %Transition Matrix is N by N
B = [.33 .33 .33 .33; .33 .33 .33 .33; .33 .33 .33 .33];%Emission Matrix is N by T

[alpha, alpham, nalpham,alphac] = alpharec(A, B, Pi);

if nalpham(2,4)/(prod(alphac)) == alpham(end,end)
    disp('succuess')
end

scalar = prod(alphac);

check = nalpham(2,4)/prod(alphac);
check2 = alpham(2,4);

[beta, betam] = betarec(A, B, Pi, alphac);

[N, T] = size(betam);

gammaa = HMMgammaposteri(1, 1, alpham, betam);

[xi,xiold] = HMMxiposteri(1, 1, 1, alpham, betam, A, B, alphac);

