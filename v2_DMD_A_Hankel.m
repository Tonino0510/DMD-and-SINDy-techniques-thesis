clear all, close all
l = 0;                                                  % Lorenz variable
params = [-1 -2 -3 -4 -5 -20 -30 -40 -50];
x0 = ones(1,9);                                         % initial condition
dt = 0.001;
tspan = 0:dt:5-dt;                                      % 5.000 time span
options = odeset('RelTol',1e-12,'AbsTol',1e-12);        % ode options
[t,x] = ode45(@(t,x)diagonal_v1(t,x,params),tspan,x0,options);

%% hankel matrix
x = x.';
x = awgn(x,10,'measured');
%plot([x y])
%legend('Original Signal','Signal with AWGN')
%x=y; %added noise
n = 20;                                                 % number of rows
m = length(tspan)-n;                                    % number of columns
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

%% DMD of x
[Phi,b,omega,r] = DMD(X,Xprime,dt,l);

%% DMD reconstruction
time_dynamics = zeros(r, m);               
for iter = 1:m
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
Xdmd = real(Phi * time_dynamics);                              
Xdmdf = Xdmd(1:n:end,:);

%% error
error=x(:,1:m)'-Xdmdf';
plot(error,'LineWidth',2)
legend('variable 1','variable 2','variable 3','variable 4','variable 5','variable 6','variable 7','variable 8','variable 9')
xlabel('sample');
ylabel('amplitude');
mse(error)
figure
plot(x(1,1:m),'-.','LineWidth',3),hold on, plot(x(2,1:m),'-.','LineWidth',3), hold on, plot(x(3,1:m),'-.','LineWidth',3),
plot(x(:,4),'-.','LineWidth',3),hold on, plot(x(:,5),'-.','LineWidth',3), hold on, plot(x(6,1:m),'-.','LineWidth',3),
plot(x(7,1:m),'-.','LineWidth',3),hold on, plot(x(8,1:m),'-.','LineWidth',3), hold on, plot(x(9,1:m),'-.','LineWidth',3),
plot(Xdmd(1,:),'--','LineWidth',3),hold on, plot(Xdmd(2,:),'--','LineWidth',3), hold on, plot(Xdmd(3,:),'--','LineWidth',3),
plot(Xdmd(4,:),'--','LineWidth',3),hold on, plot(Xdmd(5,:),'--','LineWidth',3), hold on, plot(Xdmd(6,:),'--','LineWidth',3),
plot(Xdmd(7,:),'--','LineWidth',3),hold on, plot(Xdmd(8,:),'--','LineWidth',3), hold on, plot(Xdmd(9,:),'--','LineWidth',3),
legend('x1 reale','x2 reale','x3 reale','x4 reale','x5 reale','x6 reale','x7 reale','x8 reale','x9 reale','x1 stimato','x2 stimato','x3 stimato','x4 stimato','x5 stimato','x6 stimato','x7 stimato','x8 stimato','x9 stimato')
xlabel('sample');
ylabel('amplitude')
figure 
plot(x(1,1:m),'-r','LineWidth',2), hold on, plot(Xdmdf(1,:),'--b','LineWidth',3)
legend('x1 reale','x1 stimato')
xlabel('time(s)');
ylabel('amplitude');
x_1 = x(:,1:m);
c1 = corrcoef(x_1,Xdmd_1);
c2 = c1(1,2);
