function [ wn ] = Esprit( signal, K )

    N = size(signal, 2);
    Ryy = CorrelationMatrix(signal, N-K);
    % On utilise svd car elle trie les valeurs propres.
    [U, D, V] = svd(Ryy);
    V = U(:,1:K);
    V1 = V (1:N-K-1, :);
    V2 = V (2:N-K, :);
    psi = V1\V2;
    vp = eig(psi);
    wn = angle(vp);

end
