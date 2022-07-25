function dx = diagonal_v2(t,x,A,B,g,gt)
%% Lorenz equations
g = interp1(gt,g,t);
dx = [                                         
     A(1)*x(1) + B(1)*g;
     A(2)*x(2) + B(2)*g;
     A(3)*x(3) + B(3)*g; 
     A(4)*x(4) + B(4)*g;
     A(5)*x(5) + B(5)*g;
     A(6)*x(6) + B(6)*g; 
     A(7)*x(7) + B(7)*g;
     A(8)*x(8) + B(8)*g;
     A(9)*x(9) + B(9)*g; 
];
