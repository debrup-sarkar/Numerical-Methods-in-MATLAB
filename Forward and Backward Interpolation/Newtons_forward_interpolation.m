% Lab 10
% Newton's forward difference

clc;clear;

syms x;
syms x1;
f = matlabFunction(5^x);

x = linspace(0,1,6)
y = f(x)


del_11 = y(2) - y(1);
del_12 = y(3) - y(2);
del_13 = y(4) - y(3);
del_14 = y(5) - y(4);
del_15 = y(6) - y(5);

del_21 = del_12 - del_11;
del_22 = del_13 - del_12;
del_23 = del_14 - del_13;
del_24 = del_15 - del_14;

del_31 = del_22 - del_21;
del_32 = del_23 - del_22;
del_33 = del_24 - del_23;

del_41 = del_32 - del_31;
del_42 = del_33 - del_32;

del_51 = del_42 - del_41;

x_forc = linspace(0,1,36);

h = x_forc(2) - x_forc(1);


p = (x1-x(1))/h;

y_n = y(1) + p*(del_11) + (1/factorial(2))*p*(p-1)*del_21 + ...
    (1/factorial(3))*p*(p-1)*(p-2)*(del_31) + ...
    (1/factorial(4))*p*(p-1)*(p-2)*(p-3)*(del_41) +...
    (1/factorial(5))*p*(p-1)*(p-2)*(p-3)*(p-4)*(del_51);

y_nfd = matlabFunction(y_n);
x_nfd = linspace(0,1,36);

f_interpolated = y_nfd(x_nfd);
a = f_interpolated;
%plot(x_nfd,a)
%hold on
plot(x_nfd,f(x_nfd));
hold on
plot(x_nfd,a)
