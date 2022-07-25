clear all, close all
%% Dichiarazione variabili
l = 1;                                                % Lorenz variable
Beta = [10; 28; 8/3];                                 % chaotic values
x0 = [0; 1; 20];                                      % initial condition
dt = 0.001;
tspan = 0:dt:15-dt;                                   % time span
options = odeset('RelTol',1e-12,'AbsTol',1e-12);      % ode options
[t,x] = ode45(@(t,x)lorenz(t,x,Beta),tspan,x0,options);

%% Hankel matrix
x = x.';
n = 3000;                                              % number of rows
m = length(tspan)-n;                                   % number of columns
index1 = 1:n;                                                               % from data-driven spectral analysis for c... 
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

%% DMD of x
[Phi,b,omega,r] = DMD(X,Xprime,dt,l);

%% DMD reconstruction
time_dynamics = zeros(r, m);      
for iter = 1:m                    
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
Xdmd =real(Phi * time_dynamics);         
Xdmdf = Xdmd(1:n:end,:);

%% DMD error
error = x(:,1:m)'-Xdmdf';
plot(error,'LineWidth',2)
legend('variable 1','variable 2','variable 3')
xlabel('time(s)');
ylabel('amplitude');
mse(error)
figure
plot(x(1,1:m),'-.','LineWidth',4),hold on, plot(x(2,1:m),'-.','LineWidth',4), hold on, plot(x(3,1:m),'-.','LineWidth',4),
plot(Xdmdf(1,:),':','LineWidth',4),hold on, plot(Xdmdf(2,:),':','LineWidth',4), hold on, plot(Xdmdf(3,:),':','LineWidth',4)
legend('x1 reale','x2 reale','x3 reale','x1 stimato','x2 stimato','x3 stimato')
xlabel('time(s)');
ylabel('amplitude');
c1 = corrcoef(x(:,1:m),Xdmdf)
c2 = c1(1,2);
