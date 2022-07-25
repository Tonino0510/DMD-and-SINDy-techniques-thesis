clear all, close all
load SRU_l2_l4_Dataset

Data.l4_in_n=zscore(l4_in);
Data.l4_out_n=zscore(l4_out);
%l4_in e l4_out

%% Declaration variable
dt = 0.01;                      % Time step
l = 0;
tspan = 0:dt:110;
tr = size(Data.l4_out_n,1);
tspan = tspan(1:tr);
x0 = Data.l4_out_n(1,:);

%% hankel matrix
Data.l4_out_n = Data.l4_out_n.';
n = 2000;                                             % number of rows 
m = length(tspan)-n;                                  % number of columns
index1 = 1:n;
index2 = n:n+m-1;

X = []; Xprime=[];

for ir = 1:size(Data.l4_out_n,1)
   
    % Hankel blocks ()
    c = Data.l4_out_n(ir,index1).'; r = Data.l4_out_n(ir,index2);
    H = hankel(c,r).';
    c = Data.l4_out_n(ir,index1+1).'; r = Data.l4_out_n(ir,index2+1);
    UH= hankel(c,r).';
    
    X=[X,H]; Xprime=[Xprime,UH];
end
X=X';Xprime=Xprime';

%% construction of the snapshot matrices
Ups = Data.l4_in_n(1:size(X,2),:)';
Omega = [X;Ups];
Omega1 = Omega(:,1:end-1);
OmegaPrime = Omega(:,2:end);

%% compute the svd of the input space Omega         
%for k=1:11
[Phi,b,omega,r,Ur,Sigmar,Vr] = DMD(Omega1,OmegaPrime,dt,l);

%% Recostruction
time_dynamics = zeros(r, length(tspan));
for iter = 1:length(tspan)
    time_dynamics(:,iter) = (b.*exp(omega*tspan(iter)));
end
Xdmd = real(Phi * time_dynamics);        
Xdmdf = [Xdmd(1,:)',Xdmd(n+1,:)'];

%% error
err = Data.l4_out_n'-Xdmdf;
%uscita1(k,:) = Xdmd(1,:);
%uscita2(k,:) = Xdmd(n+1,:);
figure, plot(err,'LineWidth',2)
legend('variable 1','variable 2')
xlabel('time(sample)');
ylabel('amplitude');
ms1 = mse(err);
figure 
plot(Data.l4_out_n(1,:)','-','LineWidth',3), hold on, plot(Data.l4_out_n(2,:)','-','LineWidth',3) 
plot(Xdmd(1,:)','--','LineWidth',3), hold on, plot(Xdmd(2001,:)','--','LineWidth',3)
legend('x1 reale','x2 reale','x1 stimato','x2 stimato')
xlabel('time(sample)');
ylabel('amplitude');
c1 = corrcoef(Data.l4_out_n',[Xdmd(1,:)',Xdmd(n+1,:)']);
C = c1(1,2);

%FINE

%% doppio plot
figure
subplot(1,2,1)
plot(Data.l4_out_n(1,:)','-','LineWidth',3), hold on, plot(Xdmd(1,:)','--','LineWidth',3)
xlabel('time(sample)')
ylabel('amplitude')
legend('x1 reale','x1 stimato')
subplot(1,2,2)
plot(Data.l4_out_n(2,:)','-','LineWidth',3), hold on, plot(Xdmd(2001,:)','--','LineWidth',3)
xlabel('time(sample)')
ylabel('amplitude')
legend('x2 reale','x2 stimato')

%% plot mse
%farlo sia per mse che corrcoef. citare caso migliore 50 e articolo
figure
for j=1:10
for i=1:11
    er(i,j) = mse(Data.l4_out_n(:,1:1000*j)'- [uscita1(i,1:1000*j);uscita2(i,1:1000*j)]');
end
plot(er(:,j),'LineWidth',2), hold on
end
legend('1000 sample','2000 sample','3000 sample','4000 sample','5000 sample','6000 sample','7000 sample','8000 sample','9000 sample','10000 sample')
xticklabels({'r = 5','r = 10','r = 20','r = 25','r = 50','r = 100','r = 250','r = 500','r = 1000','r = 2000','r = 3000',})
ylabel('amplitude (log)');
set(gca, 'YScale', 'log')

%% plot coefficiente di correlazione
figure
for j=1:10    
for i=1:11
    c1 = corrcoef(Data.l4_out_n(:,1:1000*j)', [uscita1(i,1:1000*j);uscita2(i,1:1000*j)]');
    cf(i,j) = c1(1,2);
end
plot(cf(:,j),'LineWidth',2), hold on
end
legend('1000 sample','2000 sample','3000 sample','4000 sample','5000 sample','6000 sample','7000 sample','8000 sample','9000 sample','10000 sample')
xticklabels({'r = 5','r = 10','r = 20','r = 25','r = 50','r = 100','r = 250','r = 500','r = 1000','r = 2000','r = 3000',})
ylabel('amplitude');