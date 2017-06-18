function [ wn ] = MatrixPencil( signal, K , p)
    N = size(signal, 2);
    Y0 = hankel(signal(1:N-p), signal(N-p:N-1));
    Y1 = hankel(signal(2:N-p+1), signal(N-p+1:N));
    [U, S, V] = svd(Y1, 0);
    [U, S, V] = sortEigenvalues(U, S, V, K);
    %Ytilde=U(1:K)*S(1:K, 1:K)*V(1:K)';
    Zp=V*inv(S)*(U')*Y0;
    Sz = eig(Zp);
    %disp(Sz(1:K+1));
    %disp(1/Sz(1:K+1));
    for k=1:K
        wn(k, 1) = -angle(1/Sz(k));
    end
    %disp(wn);
end


function [ Us, Ss, Vs ]=sortEigenvalues(U, S, V, K)
%Sort by decreasing eigenvalues 
% and take the K eigenvectors corresponding to K higher eigenvalues
    N = size(V,2);
    for i=1:N
        for j=i+1:N
            if (abs(S(i,i)) < abs(S(j,j)))
                temp = S(j,j);
                S(j,j) = S(i,i);
                S(i,i) = temp;
                temp2 = V(:,j);
                V(:,j) = V(:,i);
                V(:,i) = temp2;
                temp3 = U(:,j);
                U(:,j) = U(:,i);
                U(:,i) = temp3;
            end
        end
    end
    Us = U(:, 1:K);
    Vs = V(:, 1:K);
    Ss = S(1:K, 1:K);
end
