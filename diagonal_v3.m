function dx = diagonal_v3(t,x,A,B,g,gt)
%% Lorenz equations
g = interp1(gt,g,t);
dx = [                                         
     A(1)*x(1) + B(1)*g;
     A(2)*x(2) + B(2)*g;
     A(3)*x(3) + B(3)*g; 
     A(4)*x(4) + B(4)*g;
     A(5)*x(5) + B(5)*g;
];
