clc;clear;close;


% Gauss_Seidel
% Input matrix A and vector b
A = [4 -1 0; -1 4 -1; 0 -1 4];
b = [10; 20; 30];

% Initialize x with zeros
x = zeros(size(b));

% Set the number of iterations and tolerance
max_iter = 100;
tolerance = 1e-6;

% Gauss-Seidel iteration
for iter = 1:max_iter
    x_old = x;
    for i = 1:length(b)
        x(i) = (b(i) - A(i,1:i-1)*x(1:i-1) - A(i,i+1:end)*x_old(i+1:end)) / A(i,i);
    end
    % Check for convergence
    if norm(x - x_old) < tolerance
        break;
    end
end

% Output the result
fprintf('Solution vector: \n');
disp(x);
fprintf('Number of iterations: %d\n', iter);
