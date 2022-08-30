function [u,x,t] = heat_richardson(a,xf,T,it0,bx0,bxf,M,N)
%solve a u_xx = u_t for 0 <= x <= xf, 0 <= t <= T
% Initial Condition: u(x,0) = it0(x)
% Boundary Condition: u(0,t) = bx0(t), u(xf,t) = bxf(t)
% M = # of subintervals along x axis
% N = # of subintervals along t axis
dx = xf/M; x = [0:M]'*dx;
dt = T/N; t = [0:N]*dt;
for i = 1:M + 1, u(i,1) = it0(x(i)); end
for n = 1:N + 1, u([1 M + 1],n) = [bx0(t(n)); bxf(t(n))]; end
r = a*dt/dx/dx;

%------------------------------------------------------------
%using the explicit FDM method to guess the solution for u(i,2)
r = a*dt/dx/dx, r1 = 1 - 2*r;
for k = 1
for i = 2:M
  u(i,k+1) = r*(u(i + 1,k) + u(i-1,k)) + r1*u(i,k);
end
end


% r1 = 2*(1-r); r2 = 2*(1+r);
% for i = 1:M-1
% A(i,i) = r2; 
% if i > 1, A(i-1,i) = -r; A(i,i-1) = -r; end
% end
% for k = 2
% b = [r*u(1,k); zeros(M-3,1); r*u(M+1,k)] ...
% + r*(u(1:M-1,k-1)+u(3:M+1,k-1))+r1*u(2:M,k-1);
% u(2:M,k) = trid(A,b); 
% end

%using the explicit Richardson method to get the solution for
%u(i,3),u(i,4),u(i,5)....
for k = 1:N-1
for i = 2:M
u(i,k+2) = 2.0*r*(u(i+1,k+1) + u(i-1,k+1) -2.0*u(i,k+1)) +u(i,k); %Richardson FDM method 
end
end