function dx = diagonal_v1(t,x,A)
%% Lorenz equations
dx = [                                         
     A(1)*x(1);
     A(2)*x(2);
     A(3)*x(3); 
     A(4)*x(4);
     A(5)*x(5);
     A(6)*x(6); 
     A(7)*x(7);
     A(8)*x(8);
     A(9)*x(9); 
];
