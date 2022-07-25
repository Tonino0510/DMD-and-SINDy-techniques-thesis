clear all, close all
load SRU_l2_l4_Dataset
Data.l4_in_n=zscore(l4_in);
Data.l4_out_n=zscore(l4_out);

%% Dichiarazione variabili
dt = 0.01;                                    % Time step
l = 0;                                        % Lorenz variable
tspan = 0:dt:110;
tr = size(Data.l4_out_n,1);
tspan = tspan(1:tr);

%% construction of the snapshot matrices
X = transpose(Data.l4_out_n(1:end-1,:));
Xprime = transpose(Data.l4_out_n(2:end,:));
Ups = transpose(Data.l4_in_n(1:end-1,:));
Omega = [X;Ups];
OmegaOne = Omega(:,1:end-1);
OmegaPrime = Omega(:,2:end);

%% compute the svd of the input space Omega
[~,Ur,Sigmar,Vr] = DMDc(Omega,l);                     

%% compute the svd of the output space Xprime
[~,Uhat,Sighat,Vhat] = DMDc(Xprime,l);

%% compute the approximation of the operator G = [A B]
n = size(X,1);                  % length of the first dimension of X
q = size(Ups,1);                % length of the first dimension of Ups
U_1 = Ur(1:n,:);
U_2 = Ur(n+q:n+q,:);            

%% Evoluzione sistema discreto
approxAd = (Xprime)*Vr*inv(Sigmar)*U_1';
approxBd = (Xprime)*Vr*inv(Sigmar)*U_2';         
approxX(:,1)=Data.l4_out_n(1,:);                
for k=1:10080
approxX(:,k+1)=approxAd*approxX(:,k)+approxBd*Data.l4_in_n(k);
end

%% Error
figure 
plot(Data.l4_out_n(:,1),':','LineWidth',5),hold on, plot(Data.l4_out_n(:,2),':','LineWidth',5)
plot(approxX(1,:)','-','LineWidth',5),hold on, plot(approxX(2,:)','-','LineWidth',5)
legend('x1 reale','x2 reale','x1 stimato','x2 stimato')
xlabel('time(s)');
ylabel('amplitude');

err = Data.l4_out_n(:,:)-approxX(:,:)';
figure, plot(err,'LineWidth',2)
legend('variable 1','variable 2','variable 3')
xlabel('time(s)');
ylabel('amplitude');
ms2 = mse(err);
c1 = corrcoef(Data.l4_out_n(1:size(approxX,2),:),approxX(:,:)');
c2 = c1(1,2);
