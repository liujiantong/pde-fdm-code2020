%solve_heat
a = 1; %the parameter of (E9.2.1)
it0 = inline('sin(pi*x)','x'); %initial condition
bx0 = inline('0'); bxf = inline('0'); %boundary condition
xf = 1; M = 20; T = 0.1; N = 100; %r = 0.625
%analytical solution
uo = inline('sin(pi*x)*exp(-pi*pi*t)','x','t');
[u1,x,t] = heat_exp(a,xf,T,it0,bx0,bxf,M,N); %converge conditionally
figure(1), clf, mesh(t,x,u1)
xlabel('t');
ylabel('x');
zlabel('近似解u(x,t)');
title('最简显格式');

[u2,x,t] = heat_imp(a,xf,T,it0,bx0,bxf,M,N); %converge unconditionally
figure(2), clf, mesh(t,x,u2);
xlabel('t');
ylabel('x');
zlabel('u(x,t)');

title('最简隐格式');

[u3,x,t] = heat_CN(a,xf,T,it0,bx0,bxf,M,N); %converge unconditionally
figure(3), clf, mesh(t,x,u3);
xlabel('t');
ylabel('x');
zlabel('u(x,t)');

title('Crank-Nicolson格式');
% 
[u4,x,t] = heat_richardson(a,xf,T,it0,bx0,bxf,M,N); %does not converge 
figure(4), clf, mesh(t,x,u4);
xlabel('t');
ylabel('x');
zlabel('u(x,t)');

title('Richardson格式');

MN = M*N;
Uo = uo(x,t); aUo = abs(Uo)+eps; %values of true analytical solution
%How far from the analytical solution?
err1 = norm((u1-Uo)./aUo)/MN
err2 = norm((u2-Uo)./aUo)/MN
err3 = norm((u3-Uo)./aUo)/MN
 err4 = norm((u4-Uo)./aUo)/MN