%%%%%%%%%%% Premiere implementation %%%%%%%%%%%

% function [matCorr]=CorrelationMatrix(signal, dim)
% %Build an estimate of the correlation matrix
% % from 'signal' and of size 'dim'
%     matCorr = zeros(dim,dim);
%     for j=1:dim
%         for l=j:dim
%             rxx = correlation(signal, l-j, dim);
%             matCorr(l,j) = rxx;
%             matCorr(j,l) = conj(rxx);
%         end
%     end
% end
% 
% function [corr]=correlation(signal, ind, dim)
% %Evaluate the correlation coefficient of index 'ind'
%     N = size(signal, 2);
%     m = ind;
%     if (ind < 0)
%         m = -ind;
%     end
%     corr = signal(1:N-dim+1)*signal(1+m:N+m-dim+1)';
%     if (ind < 0)
%         corr = conj(corr);
%     end
%     corr = corr/(N);
% end

%%%%%%%%%%%%% Deuxieme implementation %%%%%%%%%%%

function [matCorr]=CorrelationMatrix(signal, dim)
    matCorr = zeros(dim,dim);
    N = size(signal,2);
    for i=0:N-dim
        matCorr = matCorr + signal(i+1:i+dim)'*signal(i+1:i+dim);
    end
    matCorr = matCorr/N;
end