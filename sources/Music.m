function r=Music(vector, p)
tol=10e-4;
%1)Calcul de la matrice de covariance
M=CorrelationMatrix(vector,n);
%2)Recherche du vecteur propre associÃ© la valeur propre minimale
[V,D]=eig(M);
s=size(D);
s=s(:,1);
Val=diag(D);
Valsort=sort(abs(Val));
%3)Calcul des polynômes annulateurs
polynome=0;
for i=1:s-p
    x=Valsort(i);
    for j=i:s
        if abs(Val(j))-x <= tol
           polynome=polynome+conv(V(:,j),conj(flipud(V(:,j))));
        end
        break;
    end
end
racines=roots(polynome);
modracines=racines(abs(racines)<1);
[tmp,idx]=sort(abs(abs(modracines)-1));
racines=modracines(idx);
freq=-angle(racines(1:p));
r=freq
end







