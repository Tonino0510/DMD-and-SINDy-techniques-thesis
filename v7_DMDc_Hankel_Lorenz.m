clear all, close all
l = 1;                                                % Lorenz variable
Beta = [8; 35; 8/3];                                  % chaotic values
x0 = [0; 1; 20];                                      % initial condition
dt = 0.001;
tspan = 0:dt:9-dt;                                    % time span
omeg = 4.5;                                           % input pulsation
a = 6;                                                % input amplitude
g = (sin(omeg*tspan));
t1=0:dt:50-dt;
options = odeset('RelTol',1e-12,'AbsTol',1e-12);      % ode options
[t,x] = ode45(@(t,x)lorenz_c(t,x,Beta,g,tspan,a),tspan,x0,options);

%% hankel matrix
x = x.';
n = 3000;                                             % number of rows  
m = length(tspan)-n;                                  % number of columns
index1 = 1:n;
index2 = n:n+m-1;

X = []; Xprime=[];

for ir = 1:size(x,1)
   
    % Hankel blocks ()
    c = x(ir,index1).'; r = x(ir,index2);
    H = hankel(c,r).';
    c = x(ir,index1+1).'; r = x(ir,index2+1);
    UH= hankel(c,r).';
    
    X=[X,H]; Xprime=[Xprime,UH];
end
X=X';Xprime=Xprime';

%% construction of the snapshot matrices
Ups = (g(:,1:size(X,2)));
Omega = [X;Ups];
Omega1 = Omega(:,1:end-1);
OmegaPrime = Omega(:,2:end);

%% compute the svd of the input space Omega

[Phi,b,omega,r,Ur,Sigmar,Vr] = DMD(Omega1,OmegaPrime,dt,l); 

%% DMD reconstruction
time_dynamics = zeros(r, m);
for iter = 1:m
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
Xdmd = real(Phi * time_dynamics);        
Xdmdf = Xdmd(1:n:end-1,:);

%% error
err = x(:,1:m)'-Xdmdf';
%uscita1(k,:) = Xdmdf(1,:);
%uscita2(k,:) = Xdmdf(2,:);
%uscita3(k,:) = Xdmdf(3,:);
figure, plot(err,'LineWidth',2)
legend('variable 1','variable 2','variable 3')
xlabel(['time(sample)',num2str(k)]);
ylabel('amplitude');
MSc= mse(err);
c1 = corrcoef(x(:,1:m)',Xdmdf');
c2 = c1(1,2);
figure 
plot(x(1,1:m),'-.','LineWidth',4),hold on, plot(x(2,1:m),'-.','LineWidth',4), hold on, plot(x(3,1:m),'-.','LineWidth',4),
plot(Xdmdf(1,:),':','LineWidth',4),hold on, plot(Xdmdf(2,:),':','LineWidth',4), hold on, plot(Xdmdf(3,:),':','LineWidth',4)
legend('x1 reale','x2 reale','x3 reale','x1 stimato','x2 stimato','x3 stimato')
xlabel('time(s)');
ylabel('amplitude');

%FINE

%% plot mse
figure
for j=1:6
for i=1:11
    er(i,j) = mse(x(:,1:1000*j)'- [uscita1(i,1:1000*j);uscita2(i,1:1000*j);uscita3(i,1:1000*j)]');
end
plot(er(:,j),'LineWidth',2), hold on
end
legend('1000 sample','2000 sample','3000 sample','4000 sample','5000 sample','6000 sample')
xticklabels({'r = 5','r = 10','r = 20','r = 25','r = 50','r = 100','r = 250','r = 500','r = 1000','r = 2000','r = 3000',})
ylabel('amplitude (log)');
set(gca, 'YScale', 'log')

%% plot coefficiente di correlazione
figure
for j=1:6    
for i=1:11
    c1 = corrcoef(x(:,1:1000*j)', [uscita1(i,1:1000*j);uscita2(i,1:1000*j);uscita3(i,1:1000*j)]');
    cf(i,j) = c1(1,2);
end
plot(cf(:,j),'LineWidth',2), hold on
end
legend('1000 sample','2000 sample','3000 sample','4000 sample','5000 sample','6000 sample')
xticklabels({'r = 5','r = 10','r = 20','r = 25','r = 50','r = 100','r = 250','r = 500','r = 1000','r = 2000','r = 3000',})
ylabel('amplitude');