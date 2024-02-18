% Simpsons one-third rule


clc;
clear;

syms x;
f = matlabFunction(1 - x - 4*(x^3) + 2*(x^5));
a = -2;
b = 4;
n = 100;

h = (b-a)/n;
x = linspace(a,b,n+1);
y = f(x);

% Use the trapezoidal rule to approximate the integral
I = h/3 * (y(1) + 4*sum(y(2:end-1) + 2*sum(y(2:end-2))) + y(end));
fprintf('Value of integration (Simpsons 1/3)\n')
disp(I)

