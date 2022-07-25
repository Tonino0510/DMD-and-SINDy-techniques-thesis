function dx = lorenz_d(t,x,Beta,g,gt,a)
%% Lorenz equations
dx = [
    -Beta(1)*(x(1)-x(2)-(a*g));
    x(1)*(Beta(2)-x(3)) - x(2);
    x(1)*x(2) - Beta(3)*x(3);
];
end