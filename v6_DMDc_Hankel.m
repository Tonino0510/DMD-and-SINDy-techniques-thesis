clear all, close all
l = 0;                                                % Lorenz variable
dt = 0.01;
tspan = 0:dt:10-dt;  
A = [-1 0 0 0 0; 0 -2 0 0 0; 0 0 -3 0 0;              % Parameter
     0 0 0 -100 0; 0 0 0 0 -100];           
B = [1 2 3 4 5]';                                     % Parameter
x0 = [1 .5 -1 2 3];                                   % Initial condition
g = transpose(3*sin(tspan-0.25));

options = odeset('RelTol',1e-12,'AbsTol',1e-12);      % ode options
[t,x] = ode45(@(t,x)diagonal_v3(t,x,A,B,g,tspan),tspan,x0,options);
%x = awgn(x,10,'measured');                           % noise

%% hankel matrix
x = x.';
n = 200                                        % number of rows 20 STAMPARE MSE E CORRCOEF PER OGNI N MODIFICATO
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
Ups = g(1:size(X,2),:)';
Omega = [X;Ups];
Omega1 = Omega(:,1:end-1);
OmegaPrime = Omega(:,2:end);

%% compute the svd of the input space Omega                   
[Phi,b,omega,r,Ur,Sigmar,Vr] = DMD(Omega1,OmegaPrime,dt,l);
%% compute the svd of the output space X'
[~,Uhat,Sighat,Vhat] = DMDc(Xprime,l); 

%% compute the approximation of the operator G = [A B]
nn = size(X,1);                  % length of the first dimension of X
q = size(Ups,1);                 % length of the first dimension of Ups
U_1 = Ur(1:nn,:);
U_2 = Ur(nn+q:nn+q,:);            

%% DMD reconstruction
time_dynamics = zeros(r, length(t));
for iter = 1:length(t)
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
Xdmd = real(Phi * time_dynamics);
Xdmdf = Xdmd(1:n:end-1,:);

%% error
err = x'-Xdmdf';
mse(err)
c1 = corrcoef(x',Xdmdf');
c2 = c1(1,2)
figure, plot(err,'LineWidth',2)
legend('variable 1','variable 2','variable 3','variable 4','variable 5')
xlabel('sample');
ylabel('amplitude');
figure
plot(x(1,:)','-','LineWidth',5),hold on, plot(x(2,:)','-','LineWidth',5), hold on, plot(x(3,:)','-','LineWidth',5),
plot(x(4,:)','-','LineWidth',5),hold on, plot(x(5,:)','-','LineWidth',5)
plot(Xdmdf(1,:)',':','LineWidth',5),hold on, plot(Xdmdf(2,:)',':','LineWidth',5), hold on, plot(Xdmdf(3,:)',':','LineWidth',5),
plot(Xdmdf(4,:)',':','LineWidth',5),hold on, plot(Xdmdf(5,:)',':','LineWidth',5)
legend('x1 reale','x2 reale','x3 reale','x4 reale','x5 reale','x1 stimato','x2 stimato','x3 stimato','x4 stimato','x5 stimato')
xlabel('sample');
ylabel('amplitude');
