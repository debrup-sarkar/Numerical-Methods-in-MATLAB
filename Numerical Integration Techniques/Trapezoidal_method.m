% Name: DEBRUP SARKAR
% Q1(a) Trapezoidal Method

clc;
clear;
close;

syms x;
f1 = 1 - x - 4*(x^3) + 2*(x^5);
f = matlabFunction(f1);

% f = @(x) ( 1 - x - 4*(x^3) + 2*(x^5)); %f(x)
a = -2; % lower limit
b = 4; % upper limit
n = 100; % number of segments

% a = input('a (lower limit) = ');
% b = input('b (Upper limit) = ');
% n = input('n (Number of segments) = ');


h = (b-a)/n;
x = linspace(a,b,n+1);
y = f(x);


I = h/2 * (y(1) + 2*sum(y(2:end-1)) + y(end));
fprintf('Value of integration (Trapezoidal Method)\n');
disp(I);




