function [r,Ur,Sigmar,Vr] = DMDc(X,lorenz)

if (lorenz==0) 
    thresh  = 1e-10;
else thresh = 1e-6;
end

[U,Sigma,V] = svd(X,'econ');                  % Step 1 con econ riduce il
r = length(find(diag(Sigma)>thresh));         % truncation value

Ur = U(:,1:r);                                % left singular vector matrix 
Sigmar = Sigma(1:r,1:r);                      % singular values matrix 
Vr = V(:,1:r);                                % right singular vector matrix 
