clc;
clear;

% Define the coefficient matrix and the right-hand side vector
A = [10 2 -1; -3 -6 2; 1 1 5];
b = [27; -61.5; -21.5];

% A = input('Enter the coeff. matrix: \n')
% b = input('Enter b(Intercept): \n')

% Performing LU decomposition of A
n = size(A, 1);
L = eye(n);
U = A;

for k = 1:n-1
    for i = k+1:n
        L(i,k) = U(i,k) / U(k,k);
        U(i,k:n) = U(i,k:n) - L(i,k) * U(k,k:n);
    end
end

y = zeros(n, 1);
x = zeros(n, 1);

% solve Ly=b (forward subs.)
for i = 1:n
    y(i) = b(i) - L(i,1:i-1)*y(1:i-1);
end

% solve Ux=y (backward subs.)
for i = n:-1:1
    x(i) = (y(i) - U(i,i+1:n)*x(i+1:n)) / U(i,i);
end

disp('The solution x is ')
disp(x)