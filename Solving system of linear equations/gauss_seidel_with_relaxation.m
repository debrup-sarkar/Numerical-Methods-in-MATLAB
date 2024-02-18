
clc;
clear;

A = [-8 1 -2 ; 2 -6 -1 ; -3 -1 7];
b = [-20 ; -38 ; -34];
tol = 0.05;
max_iter =1500;

% relaxation (lambda) = 1.2
lambda = 1.2;

%   x: solution vector
%   iter: number of iterations performed

% Check that A is a square matrix
[n, m] = size(A);
if n ~= m
    error('A must be a square matrix');
end

% Check that b has the same number of rows as A
if length(b) ~= n
    error('b must have the same number of rows as A');
end

% Initialize solution vector
x = zeros(n, 1);

% Initialize iteration counter
iter = 0;

% Perform iterations until convergence or maximum number of iterations reached
while iter < max_iter
    % Update solution vector
    for i = 1:n
        x(i) = (1 - lambda)*x(i) + lambda*(b(i) - A(i,1:i-1)*x(1:i-1) - A(i,i+1:n)*x(i+1:n))/A(i,i);
    end
    
    % Check for convergence
    r = b - A*x;
    if norm(r) < tol
        break;
    end
    
    % Update iteration counter
    iter = iter + 1;
end

% Print warning if maximum number of iterations reached
if iter == max_iter
    warning('Maximum number of iterations reached before convergence');
end
fprintf('Solution(x):\n')
disp(x)
fprintf('iterations: %d \n',iter)
