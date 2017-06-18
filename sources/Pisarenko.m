
function r=Pisarenko(vector, p)
tol=10^-7;
%1)Calcul de la matrice de covariance
M=CorrelationMatrix(vector,p+1);
%2)Recherche du vecteur propre associé la valeur propre minimale
[V,D]=eig(M);
s=size(D);
s=s(:,1);
minval=D(1,1);
multiplicite=1;
for i=2:s
    if abs(D(i,i))<abs(minval)
        minval=D(i,i);
        multiplicite=1;
    elseif abs(D(i,i)-minval)<=tol
        multiplicite=multiplicite+1;
    end
end
%3)Modificaton de la matrice en fonction de la multiplicité
if multiplicite==1
    multiplicite=0;
end
M=M-minval*eye(p+1,p+1);
M=M(1:p+1-multiplicite,1:p+1-multiplicite);
%4) Recherche de la valeur propre associée à 0
[V,D]=eig(M);
s=size(V);
zerovect=V(:,1);
s=s(:,1);
for i=2:s
    if abs(D(i,i)) <=tol %Problème d'arrondi numerique
        zerovect=V(:,i);
    end
end
%5)Calcul des racines du polynome associé au vecteur propres=size(zerovect);
s=s(:,1);
polynome=flipud(zerovect);
racines=roots(polynome);
r=-angle(racines);
end
