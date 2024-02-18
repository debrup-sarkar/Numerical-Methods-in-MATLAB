clc;
clear;
close;

% Define a positive definite matrix
A = [4 12 -16; 12 37 -43; -16 -43 98]

% Perform Cholesky decomposition
n = size(A,1);
L = zeros(n,n);
for j = 1:n
    L(j,j) = sqrt(A(j,j) - sum(L(j,1:j-1).^2));
    for i = j+1:n
        L(i,j) = (A(i,j) - sum(L(i,1:j-1).*L(j,1:j-1))) / L(j,j);
    end
end

% Print the results
disp('The Cholesky decomposition of A is:');
disp(L);
