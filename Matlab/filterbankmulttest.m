A = [1 1 1 1 1 1; 2 2 2 2 2 2; 3 3 3 3 3 3;4 4 4 4 4 4; 5 5 5 5 5 5; 6 6 6 6 6 6; 7 7 7 7 7 7];
B = [1 2 3 4 5 6; 1 2 3 4 5 6; 1 2 3 4 5 6;1 2 3 4 5 6; 1 2 3 4 5 6; 1 2 3 4 5 6; 1 2 3 4 5 6];
[arows,acols] = size(A);
[brows,bcols] = size(B);
C = zeros(acols, bcols);

for P_idx = 1:acols
    for mel_idx = 1:bcols
        C(mel_idx,P_idx) = sum(A(:,P_idx).*B(:,mel_idx));
    end
end

C = log(C);

mfcc = zeros(1, acols);

for P_idx = 1:acols
   % Take the DCT of the log filterbank energies.
   C(:,P_idx) =  dct(C(:,P_idx));
   % Keep DCT coefficients 2-13, discard the rest.
   mfcc(1,P_idx) = C(1, P_idx);
end

% calculate delta features with linear least squared first-order estimator
deltas = zeros(1, acols);
M = 2;
syms t
for idx = 3:acols - 2
    for j = 1:bcols
        deltas(j, idx) = symsum(t*mfcc(j,idx+t), -M, M)/symsum(t^2, -M, M);
    end
end
