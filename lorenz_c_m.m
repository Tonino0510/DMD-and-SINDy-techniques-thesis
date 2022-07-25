function dx = lorenz_c_m(t,x,Beta,g,h,gt,a,b)
%% Lorenz equations
g = interp1(gt,g,t);
h = interp1(gt,h,t);
dx = [
    -Beta(1)*(x(1)-x(2)-(a*g));
    x(1)*(Beta(2)-x(3)) - x(2);
    x(1)*x(2) - Beta(3)*x(3) + b*h;
];
end