function [Phi,b,omega,r,Ur,Sigmar,Vr] = DMD(X,Xprime,dt,lorenz)

if (lorenz==0) 
    thresh  = 1e-10;
else thresh = 1e-6;
end

[U,Sigma,V] = svd(X,'econ');                  % Step 1 con econ riduce il
                                              % rango della matrice
r = length(find(diag(Sigma)>thresh));         % truncation value

Ur = U(:,1:r);                                % left singular vector matrix 
Sigmar = Sigma(1:r,1:r);                      % singular values matrix 
Vr = V(:,1:r);                                % right singular vector
                                              % matrix 

Atilde = Ur'*Xprime*Vr/Sigmar;                % Step 2. It's the best-fit
                                              % linear model
[W,Lambda] = eig(Atilde);                     % Step 3. They are the 
                                              % eigenvectors and the
                                              % eigenvalues
Phi = Xprime*(Vr/Sigmar)*W;   

lambda = diag(Lambda);                        % discrete-time eigenvalues
omega = log(lambda)/dt;                       % continuous-time eigenvalues

x1 = X(:, 1);       
b = pinv(Phi)*x1;                             % initial amplitude of each
                                              % mode vector