clear all
clc

g=9.81;m=68.1;cd=0.25;


h=input('Enter the step size, h: ');
t_i=input('Enter the initial time: ');
t_f=input('Enter the final time: ');
x_i=input('Enter the initial displacement: ');
v_i=input('Enter the initial velocity: ');
t= (t_i:h:t_f)';
n=length(t);

x_euler=zeros(n,1);x_heun=zeros(n,1);x_heun_with_iteration=zeros(n,1);x_RK2=zeros(n,1);x_RK4=zeros(n,1);
x_euler(1)=x_i;x_heun(1)=x_i;x_heun_with_iteration(1)=x_i;x_RK2(1)=x_i;x_RK4(1)=x_i;
v_euler=zeros(n,1);v_heun=zeros(n,1);v_heun_with_iteration=zeros(n,1);v_RK2=zeros(n,1);v_RK4=zeros(n,1);
v_euler(1)=v_i;v_heun(1)=v_i;v_heun_with_iteration(1)=v_i;v_RK2(1)=v_i;v_RK4(1)=v_i;
a_euler=zeros(n,1);a_heun=zeros(n,1);a_heun_with_iteration=zeros(n,1);a_RK2=zeros(n,1);a_RK4=zeros(n,1);
a_euler(1)= g-((cd*v_euler(1)^2)/m);
a_heun(1)= g-((cd*v_heun(1)^2)/m);
a_heun_with_iteration(1)= g-((cd*v_heun_with_iteration(1)^2)/m);
a_RK2(1)= g-((cd*v_RK2(1)^2)/m);
a_RK4(1)= g-((cd*v_RK4(1)^2)/m);

difrnc_x=zeros(n,1);difrnc_v=zeros(n,1);difrnc_a=zeros(n,1);
iter=1;epsilon=0.01;maxiter=100;difrnc_x(1)=0;difrnc_v(1)=0;

x_exact=((m/cd)*(log(cosh(t.*sqrt(g*cd/m)))));
v_exact=((sqrt(g*m/cd))*(tanh(t.*sqrt(g*cd/m))));


for i=1:n-1
    x_euler(i+1)=x_euler(i)+v_euler(i)*h;
    v_euler(i+1)=v_euler(i)+a_euler(i)*h;
    a_euler(i+1)= g-(cd*v_euler(i+1)^2)/m;
end

for i=1:n-1
    x_p_heun=x_heun(i)+v_heun(i)*h;
    v_p_heun=v_heun(i)+a_heun(i)*h;
    a_p_heun= g-(cd*v_p_heun^2)/m;
    x_heun(i+1)=x_heun(i)+(v_heun(i)+v_p_heun)*(h/2);
    v_heun(i+1)=v_heun(i)+(a_heun(i)+a_p_heun)*(h/2);
    a_heun(i+1)= g-(cd*v_heun(i+1)^2)/m;
end

for i=1:n-1
    x_p_heun_with_iteration=x_heun_with_iteration(i)+v_heun_with_iteration(i)*h;
    v_p_heun_with_iteration=v_heun_with_iteration(i)+a_heun_with_iteration(i)*h;
    a_p_heun_with_iteration= g-(cd*v_p_heun_with_iteration^2)/m;
    x_heun_with_iteration(i+1)=x_heun_with_iteration(i)+(v_heun_with_iteration(i)+v_p_heun_with_iteration)*(h/2);
    v_heun_with_iteration(i+1)=v_heun_with_iteration(i)+(a_heun_with_iteration(i)+a_p_heun_with_iteration)*(h/2);
    a_heun_with_iteration(i+1)= g-(cd*v_heun_with_iteration(i+1)^2)/m;
    difrnc_x(i+1)=(abs(x_heun_with_iteration(i+1)-x_p_heun_with_iteration)/x_heun_with_iteration(i+1))*100;
    difrnc_v(i+1)=(abs(v_heun_with_iteration(i+1)-v_p_heun_with_iteration)/v_heun_with_iteration(i+1))*100;
    while (difrnc_x(i+1)>=epsilon && difrnc_v(i+1)>=epsilon && iter<=maxiter)
        x_heun_with_iteration(i+1)=x_heun_with_iteration(i)+(v_heun_with_iteration(i)+v_p_heun_with_iteration)*(h/2);
        v_heun_with_iteration(i+1)=v_heun_with_iteration(i)+(a_heun_with_iteration(i)+a_p_heun_with_iteration)*(h/2);
        a_heun_with_iteration(i+1)= g-(cd*v_heun_with_iteration(i+1)^2)/m;
        difrnc_x(i+1)=(abs(x_heun_with_iteration(i+1)-x_p_heun_with_iteration)/x_heun_with_iteration(i+1))*100;
        difrnc_v(i+1)=(abs(v_heun_with_iteration(i+1)-v_p_heun_with_iteration)/v_heun_with_iteration(i+1))*100;
        iter=iter+1;
    end
end

for i=1:n-1
    k1_x_RK2=v_RK2(i);
    k1_v_RK2=a_RK2(i);
    x_k1_RK2=x_RK2(i)+k1_x_RK2*h;
    v_k1_RK2=v_RK2(i)+k1_v_RK2*h;
    a_k1_RK2= g-(cd*v_k1_RK2^2)/m;
    k2_x_RK2=v_k1_RK2;
    k2_v_RK2=a_k1_RK2;
    x_RK2(i+1)=x_RK2(i)+(k1_x_RK2+k2_x_RK2)*(h/2);
    v_RK2(i+1)=v_RK2(i)+(k1_v_RK2+k2_v_RK2)*(h/2);
end

for i=1:n-1
    k1_x=v_RK4(i);
    k1_v=a_RK4(i);
    x_k1=x_RK4(i)+k1_x*(h/2);
    v_k1=v_RK4(i)+k1_v*(h/2);
    a_k1= g-(cd*v_k1^2)/m;
    k2_x=v_k1;
    k2_v=a_k1;
    x_k2=x_RK4(i)+k2_x*(h/2);
    v_k2=v_RK4(i)+k2_v*(h/2);
    a_k2= g-(cd*v_k2^2)/m;
    k3_x=v_k2;
    k3_v=a_k2;
    x_k3=x_RK4(i)+k3_x*h;
    v_k3=v_RK4(i)+k3_v*h;
    a_k3= g-(cd*v_k3^2)/m;
    k4_x=v_k3;
    k4_v=a_k3;
    x_RK4(i+1)=x_RK4(i)+(k1_x+2*(k2_x+k3_x)+k4_x)*(h/6);
    v_RK4(i+1)=v_RK4(i)+(k1_v+2*(k2_v+k3_v)+k4_v)*(h/6);
end

result_x=[x_exact x_euler x_RK4];
result_v=[v_exact v_euler v_RK4];
result_x
result_v

subplot(2,1,1)
plot(t,result_x)
grid on
title ('ODE METHOD')
xlabel('Time, t (sec)')
ylabel ('Displacement, x (m)')
legend('Exact','Eulear','RK4thorder','Location','NorthEastOutside')
subplot(2,1,2)
plot(t,result_v)
grid on
title ('ODE METHOD')
xlabel('Time, t (sec)')
ylabel ('Velocity, v (m/sec)')
legend('Exact','Eulear','RK4thorder','Location','NorthEastOutside')