
A = [15 -3 -1; -3 18 -6; -4 -1 12];
b = [3800; 1200; 2350];
x0 = zeros(length(b),1);
tol = 0.05;
maxiter = 500;

% GAUSS_JACOBI Solves linear system Ax = b using the Gauss-Jacobi method.
% Inputs:
% A - coefficient matrix
% b - right-hand side vector
% x0 - initial guess for solution vector
% tol - tolerance for convergence
% maxiter - maximum number of iterations
% Outputs:
% x - solution vector
% iter - number of iterations performed

n = length(b);
x = x0;
iter = 0;
converged = false;

while ~converged && iter < maxiter
    x_old = x;
    for i = 1:n
        sigma = 0;
        for j = 1:n
            if j ~= i
                sigma = sigma + A(i,j) * x_old(j);
            end
        end
        x(i) = (b(i) - sigma) / A(i,i);
    end
    iter = iter + 1;
    converged = norm(x - x_old) < tol;
end

if ~converged
    warning('Gauss-Jacobi did not converge');
end

disp('Solution\n')
x
disp('iterations')
iter
