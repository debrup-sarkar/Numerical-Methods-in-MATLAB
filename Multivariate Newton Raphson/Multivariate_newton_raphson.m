clc;
clear;

% Multivariate Newton Raphson Algorithm for Solving system of equations


% Initialize x and y

x0 = 1.5;
y0 = 3.5;


% System of equations

% -------------------------------------------
u = x0^2 + x0*y0 - 10;
v = y0 + 3*x0*(y0^2) - 57;
% -------------------------------------------

% Form the Jacobian Matrix

% --------------------------------------------
du_dx = 2*x0 + y0;
du_dy = x0;
dv_dx = 3*(y0^2);
dv_dy = 1 + 6*x0*y0;

J = du_dx*dv_dy - du_dy*dv_dx;

% --------------------------------------------


e = 1;
tol = 0.01;

c = 0;

while e>tol
    x1 = x0 - ((u*dv_dy - v*du_dy)/J);
    y1 = y0 - ((v*du_dx - u*dv_dx)/J);
    e = abs((x1-x0)/x1)*100;
    x0 = x1;
    y0 = y1;

    u = x0^2 + x0*y0 - 10;
    v = y0 + 3*x0*(y0^2) - 57;

    du_dx = 2*x0 + y0;
    du_dy = x0;
    dv_dx = 3*(y0^2);
    dv_dy = 1 + 6*x0*y0;

    J = du_dx*dv_dy - du_dy*dv_dx;

    c = c+1;
end
fprintf('x = %f ; y = %f ; Number of iterations = %d',x0,y0,c);