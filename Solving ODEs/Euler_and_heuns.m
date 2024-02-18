clc;
clear all;
close all;

% Define the ODE
f = @(t,y) 4*exp(0.8*t) - 0.5*y;

% Define the initial conditions
t0 = 0;
y0 = 2;
t = [0,4];

% Define the step size and the number of steps
h = 0.25;
N = (t(2)-t(1))/h;

% Initialize the arrays for t, y, and the error
t = zeros(1,N+1);
y_euler = zeros(1,N+1);
y_heun = zeros(1,N+1);
y_exact = zeros(1,N+1);

% Set the initial values
t(1) = t0;
y_euler(1) = y0;
y_heun(1) = y0;
y_exact(1) = y0;

% Compute the numerical solution using Euler's method
for i=1:N
    y_euler(i+1) = y_euler(i) + h*f(t(i),y_euler(i));
    t(i+1) = t(i) + h;
end

% Compute the numerical solution using Heun's method
for i=1:N
    % Predictor step
    y_pred = y_heun(i) + h*f(t(i),y_heun(i));
    % Corrector step
    y_heun(i+1) = y_heun(i) + h/2*(f(t(i),y_heun(i)) + f(t(i+1),y_pred));
end

% Compute the analytical solution
for i=1:N+1
    y_exact(i) =  (4/1.3)*(exp(0.8*t(i)) - exp(-0.5*t(i))) + 2*exp(-0.5*t(i));
end

% Compute the error
error_euler = abs(y_exact - y_euler);
error_heun = abs(y_exact - y_heun);

% Plot the numerical solution and the error
plot(t,y_euler,'r-o',t,y_heun,'b-s',t,y_exact,'g--');
legend('Euler''s method','Heun''s method','Analytical solution');
xlabel('t');
ylabel('y');
title('Numerical solutions for ODE');
grid on;

figure;
plot(t,error_euler,'r-o',t,error_heun,'b-s');
legend('Euler''s method','Heun''s method');
xlabel('t');
ylabel('Error');
title('Errors in numerical solutions for ODE');
grid on;
