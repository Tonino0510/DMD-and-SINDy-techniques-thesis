clear all, close all
l = 0; 
params = [-1 -2 -3 -4 -5 -20 -30 -40 -50 ];
x0 = ones(1,9);                                     % initial condition
dt = 0.001;
tspan = 0:dt:10-dt;                                 % time span
options = odeset('RelTol',1e-12,'AbsTol',1e-12);    % ode options
[t,x] = ode45(@(t,x)diagonal_v1(t,x,params),tspan,x0,options);
%% Matrix creation 
X = transpose(x(1:end-1,:));
Xprime = transpose(x(2:end,:));

%% DMD
[Phi,b,omega,r] = DMD(X,Xprime,dt,l);

%% DMD reconstruction
time_dynamics = zeros(r, length(t));
for iter = 1:length(t)
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
Xdmd = real(Phi * time_dynamics);                              

%% error
error=x-Xdmd';
plot(error,'LineWidth',2)
legend('variable 1','variable 2','variable 3','variable 4','variable 5','variable 6','variable 7','variable 8','variable 9');
ylabel('amplitude') ;
xlabel('sample');
mse(error)
figure
plot(x(:,1),'-.','LineWidth',3),hold on, plot(x(:,2),'-.','LineWidth',3), hold on, plot(x(:,3),'-.','LineWidth',3),
plot(x(:,4),'-.','LineWidth',3),hold on, plot(x(:,5),'-.','LineWidth',3), hold on, plot(x(:,6),'-.','LineWidth',3),
plot(x(:,7),'-.','LineWidth',3),hold on, plot(x(:,8),'-.','LineWidth',3), hold on, plot(x(:,9),'-.','LineWidth',3),
plot(Xdmd(1,:)','--','LineWidth',3),hold on, plot(Xdmd(2,:)','--','LineWidth',3), hold on, plot(Xdmd(3,:)','--','LineWidth',3),
plot(Xdmd(4,:)','--','LineWidth',3),hold on, plot(Xdmd(5,:)','--','LineWidth',3), hold on, plot(Xdmd(6,:)','--','LineWidth',3),
plot(Xdmd(7,:)','--','LineWidth',3),hold on, plot(Xdmd(8,:)','--','LineWidth',3), hold on, plot(Xdmd(9,:)','--','LineWidth',3),
legend('x1 reale','x2 reale','x3 reale','x4 reale','x5 reale','x6 reale','x7 reale','x8 reale','x9 reale','x1 stimato','x2 stimato','x3 stimato','x4 stimato','x5 stimato','x6 stimato','x7 stimato','x8 stimato','x9 stimato')
xlabel('sample');
ylabel('amplitude');

c1 = corrcoef(x,Xdmd');
c2 = c1(1,2);
