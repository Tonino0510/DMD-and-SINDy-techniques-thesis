clear all, close all

l = 0;                                                % Lorenz variable
dt = 0.01;
tspan = 0:dt:10-dt;  
A = [-1 0 0 0 0; 0 -2 0 0 0; 0 0 -3 0 0;
     0 50 0 -100 0; 0 0 10 0 -100];
B = [1 2 3 4 5]';
x0 = [1 .5 -1 2 3];                                    % initial condition                           
g = transpose(15*cos(5*tspan-0.25));                   % control

options = odeset('RelTol',1e-12,'AbsTol',1e-12);       % ode options
[t,x] = ode45(@(t,x)System_v2(t,x,A,B,g,tspan),tspan,x0,options);

%% construction of the snapshot matrices
X = transpose(x(1:end-1,:));
Xprime = transpose(x(2:end,:));
Ups = transpose(g(1:end-1,:));
Omega = [X;Ups];
OmegaOne = Omega(:,1:end-1);
OmegaPrime = Omega(:,2:end);

%% compute the svd of the input space Omega                     
[Phi,b,omega,r,Ur,Sigmar,Vr] = DMD(OmegaOne,OmegaPrime,dt,l);

%% compute the svd of the output space X'
[~,Uhat,Sighat,Vhat] = DMDc(Xprime,l);

%% compute the approximation of the operator G = [A B]
n = size(X,1);                  % length of the first dimension of X
q = size(Ups,1);                % length of the first dimension of Ups
U_1 = Ur(1:n,:);
U_2 = Ur(n+q:n+q,:);            

%% DMD reconstruction
time_dynamics = zeros(r, length(t));
for iter = 1:length(t)
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
Xdmd = real(Phi * time_dynamics); 
Xdmdf = Xdmd(1:5,:);

err = x-Xdmdf';
figure, plot(err,'LineWidth',2)
legend('variable 1','variable 2','variable 3','variable 4','variable 5')
xlabel('sample');
ylabel('amplitude');
mse(err)
figure
plot(x(:,1),'-','LineWidth',5),hold on, plot(x(:,2),'-','LineWidth',5), hold on, plot(x(:,3),'-','LineWidth',5),
plot(x(:,4),'-','LineWidth',5),hold on, plot(x(:,5),'-','LineWidth',5),
plot(Xdmdf(1,:)',':','LineWidth',5),hold on, plot(Xdmdf(2,:)',':','LineWidth',5), hold on, plot(Xdmdf(3,:)',':','LineWidth',5),
plot(Xdmdf(4,:)',':','LineWidth',5),hold on, plot(Xdmdf(5,:)',':','LineWidth',5)
legend('x1 reale','x2 reale','x3 reale','x4 reale','x5 reale','x1 stimato','x2 stimato','x3 stimato','x4 stimato','x5 stimato')
xlabel('sample');
ylabel('amplitude');
