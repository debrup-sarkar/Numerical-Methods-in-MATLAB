A = [2 8 10; 8 4 5; 10 5 7];
A = inv(A);
x_init = ones(3,1);

tol = 0.0001;

c = 0;
e = inf;
while e>tol
    x1 = A*x_init;
    l = max(x1);
    x1 = x1/l;
    e = (norm(x1-x_init))/(norm(x1));
    x_init = x1;
    c = c +1;
end
disp('Smallest Eigen Value')
l