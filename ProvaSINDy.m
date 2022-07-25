clear all, close all, clc
%% Generate Data
usesine = 0;                                % sin knob
Beta = [10; 28; 8/3];                       % Lorenz's parameters (chaotic)
n = 3;                                      % variables number
polyorder = 4;
x0=[-8; 8; 27];                             % Initial condition  
dt = 0.001;
tspan=0:dt:15-dt;
omeg = 4.5;
a = 6;
g = (sin(omeg*tspan));
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));  
[t,x1] = ode45(@(t,x)lorenz_c(t,x,Beta,g,tspan,a),tspan,x0,options);

%% Compute true derivative
eps = 0.0;                                          % Noise level
for i=1:length(x1)
    dx(i,:) = lorenz_d(0,x1(i,:),Beta,g(i),tspan,a);   % Run Lorenz model to
                                                    % compute derivatives
end

dx(:,n+1) = 0*dx(:,n);                              % Add additional column
                                                    % of zeros for 
                                                    % derivatives of u
dx = dx + eps*randn(size(dx));                      % Add measurement noise
x1 = [x1 g'];                                       % Stack state and input
                                                    % measurements

%% Build library and compute sparse regression
Theta = poolData(x1,(n+1),polyorder,usesine);       % Generate library
lambda = 1e-1;                                      % sparsification 
                                                    % hyperparameter
Xi = sparsifyDynamics(Theta,dx,lambda,(n+1)); 
poolDataLIST({'x','y','z','u'},Xi,n+1,polyorder,usesine);

%% SINDy reconstruction
p.ahat = Xi(:,1:n);    
p.polyorder = polyorder; 
p.usesine = usesine;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[~,xMT]=ode45(@(t,x)sparseGalerkinControl(t,x,g,tspan,p),tspan,x0,options);


%% error
err = x1(:,1:end-1)-xMT;
figure
plot(err,'LineWidth',2)
legend('variable 1','variable 2','variable 3')
xlabel('sample');
ylabel('amplitude');
mse(err)
plot(x1(:,1),'-','LineWidth',2),hold on, plot(x1(:,2),'-','LineWidth',2), hold on, plot(x1(:,3),'-','LineWidth',2),
plot(xMT(:,1),':','LineWidth',6),hold on, plot(xMT(:,2),':','LineWidth',6), hold on, plot(xMT(:,3),':','LineWidth',6)
legend('x1 reale','x2 reale','x3 reale','x1 stimato','x2 stimato','x3 stimato')
xlabel('sample');
ylabel('amplitude');

%% errore2
plot(x1(:,1),'-','LineWidth',2),hold on, plot(x1(:,2),'-','LineWidth',2), hold on, plot(x1(:,3),'-','LineWidth',2),
plot(xMT(:,1),':','LineWidth',6),hold on, plot(xMT(:,2),':','LineWidth',6), hold on, plot(xMT(:,3),':','LineWidth',6)
legend('x1 reale','x2 reale','x3 reale','x1 stimato','x2 stimato','x3 stimato')
xlabel('time(s)');
ylabel('amplitude');

%% plot
figure
subplot(1,2,1)
plot3(xMT(:,1),xMT(:,2),xMT(:,3))
view(27,16)
grid on
xlabel('x','FontSize',13)
ylabel('y','FontSize',13)
zlabel('z','FontSize',13)
set(gca,'FontSize',13)
subplot(1,2,2)
plot3(x1(:,1),x1(:,2),x1(:,3))
view(27,16)
grid on
xlabel('x','FontSize',13)
ylabel('y','FontSize',13)
zlabel('z','FontSize',13)
set(gca,'FontSize',13)