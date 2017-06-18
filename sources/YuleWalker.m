function [ wn ] = YuleWalker( signal, K )
    % Implementation of Yule-Walker without estimation.
    alphas = correlationMat(signal, K)\correlationVec(signal, K);
    coefs = zeros(K+1,1);
    coefs(K+1)=1;
    for n=1:K
        coefs(n) = alphas(K+1-n);
    end
    racines = roots(coefs);
    wn = angle(racines);
end

function [corr]=correlation(signal, K, ind1, ind2)
    N = size(signal, 2);
    corr = signal(1+K-ind1:N-ind1)*signal(1+K-ind2:N-ind2)';
end

function [matCorr]=correlationMat(signal, K)
    matCorr = zeros(K);
    for j=1:K
        for l=j:K
            rxx_jl = correlation(signal, K, j, l);
            matCorr(l,j) = rxx_jl;
            matCorr(j,l) = conj(rxx_jl);
        end
    end
end

function [vecCorr]=correlationVec(signal, K)
    vecCorr = zeros(K,1);
    for n=1:K
        vecCorr(n) = -correlation(signal, K, 0, n);
    end
end
