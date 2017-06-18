% Prony TLS method extracted from "pulses-recovery.m"
% toolbox by Laurent Condat.

function [wn] = PronyTLS (input, K)
    N = size(input, 2);
    Tnoisy=toeplitz(input(K+1:N),input(K+1:-1:1));
    [U S V]=svd(Tnoisy,0);
    wn=-angle(roots(V(:,K+1)));	