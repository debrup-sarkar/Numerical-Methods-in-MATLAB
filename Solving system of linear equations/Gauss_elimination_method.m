clc;
clear;
close;

%C = input('Enter the coefficient matrix(C): ');
%b = input('Enter vector b :' );

C = [3 -0.1 -0.2 1 1.2; 0.1 7 -0.3 2 2.1; 0.3 -0.2 10 5 2.2; 2 1 4 2 1; 1 2 1 3.2 1.3];
b= [7.85 -19.3 71.4 1.2 10]';
A = [C b];   %Augmented Matrix
n= size(A,1);   %number of equations
x = zeros(n,1); 
for i=1:n-1
    for j=i+1:n
        m = A(j,i)/A(i,i);
        A(j,:) = A(j,:) - m*A(i,:);
    end
end
x(n) = A(n,n+1)/A(n,n);
for i=n-1:-1:1
    summ = 0;
    for j=i+1:n
        summ = summ + A(i,j)*x(j,:);
        x(i,:) = (A(i,n+1) - summ)/A(i,i);
    end
end
C
b

x
x_anal = C\b
% A