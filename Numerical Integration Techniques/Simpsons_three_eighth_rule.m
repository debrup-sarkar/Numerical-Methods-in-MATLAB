clc;
clear;

f = @(x) 1 - x - 4*(x^3) + 2*(x^5); % function to integrate
a = -2; % lower limit
b = 4; % upper limit
n = 120; % number of segments
I = simpson38(f, a, b,n); % approximate integral using Simpson's 3/8 rule
disp(I); % display the result

function I = simpson38(f, a, b,n)
% Simpon's 3/8 rule to approximate the integral of f from a to b

h = (b-a) / n; % step size
x1 = a + h;
x2 = a + 2*h;

I = (b-a) / 8 * (f(a) + 3*f(x1) + 3*f(x2) + f(b)); % 3/8 rule formula
end
