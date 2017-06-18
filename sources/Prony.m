function [wn] = Prony (input, K)
    N = size(input, 2);
    DiffMat=toeplitz(input(K:N-1), fliplr(input(1:K)));
    DiffRes=input(K+1:N);
    hj=-linsolve(DiffMat, DiffRes.');
    hj =[1; hj(1:K)];
    wn = -angle(roots(hj));