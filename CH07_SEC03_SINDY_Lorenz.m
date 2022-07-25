clear all, close all, clc
%% Generate Data
usesine = 0;                            % Using sine functions in library
Beta = [10; 28; 8/3];                   % Lorenz's parameters (chaotic)
n = 3;                                  % variables number
x0=[-8; 8; 27];                         % Initial condition
dt = 0.01;
tspan=[0.01:dt:50];
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));  
[t,x]=ode45(@(t,x) lorenz(t,x,Beta),tspan,x0,options);

%% Compute Derivative
for i=1:length(x)
    dx(i,:) = lorenz(0,x(i,:),Beta);
end

%% Build library and compute sparse regression
Theta = poolData(x,n,4,usesine);                    % Generate library
lambda = 0.025;                                     % sparsification knob
Xi = sparsifyDynamics(Theta,dx,lambda,n);
poolDataLIST({'x','y','z'},Xi,n,4);

%% Reconstruction for T\in[0,20]
tspan = [0 20];
[tA,xA]=ode45(@(t,x)lorenz(t,x,Beta),tspan,x0,options);                  % true model
[tB,xB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,4,usesine),tspan,x0,options);  % approximate

figure
subplot(1,2,1)
dtA = [0; diff(tA)];
color_line3(xA(:,1),xA(:,2),xA(:,3),dtA,'LineWidth',2);
view(27,16)
grid on
xlabel('x','FontSize',13)
ylabel('y','FontSize',13)
zlabel('z','FontSize',13)
set(gca,'FontSize',13)
subplot(1,2,2)
dtB = [0; diff(tB)];
color_line3(xB(:,1),xB(:,2),xB(:,3),dtB,'LineWidth',2);
view(27,16)
grid on
xlabel('x','FontSize',13)
ylabel('y','FontSize',13)
zlabel('z','FontSize',13)
set(gca,'FontSize',13)

%% Lorenz for t=50, dynamo view
figure
subplot(1,2,1)
plot(tA,xA(:,1),'k','LineWidth',1.5), hold on
plot(tB,xB(:,1),'r--','LineWidth',1.5)
grid on
xlabel('Time','FontSize',13)
ylabel('x','FontSize',13)
set(gca,'FontSize',13)
subplot(1,2,2)
plot(tA,xA(:,2),'k','LineWidth',1.5), hold on
plot(tB,xB(:,2),'r--','LineWidth',1.5)
grid on

%% compute error
err=xA-xB;
figure, plot(err,'LineWidth',2)
legend('variable 1','variable 2','variable 3')
xlabel('time(s)');
ylabel('amplitude');
mse(err)
figure
plot(xA(:,1),'-','LineWidth',6),hold on, plot(xA(:,2),'-','LineWidth',6), hold on, plot(xA(:,3),'-','LineWidth',6),
plot(xB(:,1),':','LineWidth',6),hold on, plot(xB(:,2),':','LineWidth',6), hold on, plot(xB(:,3),':','LineWidth',6)
legend('x1 reale','x2 reale','x3 reale','x1 stimato','x2 stimato','x3 stimato')
xlabel('time(s)');
ylabel('amplitude');

%% FIGURE 2: LORENZ for T\in[0,250]
tspan = [0 250];
options = odeset('RelTol',1e-6,'AbsTol',1e-6*ones(1,n));
[tA,xA]=ode45(@(t,x)lorenz(t,x,Beta),tspan,x0,options);   % true model
[tB,xB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,4,usesine),tspan,x0,options);  % approximate

figure
subplot(1,2,1)
dtA = [0; diff(tA)];
color_line3(xA(:,1),xA(:,2),xA(:,3),dtA,'LineWidth',1.5);
view(27,16)
grid on
xlabel('x','FontSize',13)
ylabel('y','FontSize',13)
zlabel('z','FontSize',13)

subplot(1,2,2)
dtB = [0; diff(tB)];
color_line3(xB(:,1),xB(:,2),xB(:,3),dtB,'LineWidth',1.5);
view(27,16)
grid on
xlabel('x')
ylabel('y')
zlabel('z')

