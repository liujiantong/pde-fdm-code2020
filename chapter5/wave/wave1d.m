function wave1d(N,M,T,sigma)
%
% Solve the wave equation on a 1D domain with zero dirichlet BC
% using leapfrog
%
% N = Number of spatial subintervals so h=1/N
% M = Number of time steps per unit time so dt=1/M
% T = final time
% sigma=conductivity
dt=1/M;
h=1/N;
disp(sprintf('CFL parameter dt/h = %g',dt/h))
%
% Finite difference in time matrix
Ah=makematrix(N);
% Initial condition is u=ix, u_t=0
x=linspace(0,1,N+1)';
u0=ix(x);
u0(1)=0;
u0(N+1)=0;
u=u0;
uprev=u0;
plot(x,u);%axis([0,1,-1.5,1.5])

fig=figure;
       set(fig,'DoubleBuffer','on');
       set(gca,...
       	   'NextPlot','replace','Visible','off')
       mov = avifile('example.avi')
       
Mf(1)=getframe;
unew=zeros(N+1,1);
% Courant parameter
lambda=dt^2/h^2;
for j=1:ceil(M*T)
unew(2:N)=(2*u(2:N)-(1-dt*sigma/2)*uprev(2:N)-lambda*Ah*u(2:N))/(1+sigma*dt/2);
uprev=u;
u=unew;
h=plot(x,u);%axis([0,1,-1.5,1.5]);
set(h,'EraseMode','xor');
F = getframe(gca);
%Mf(j+1)=getframe;
F = getframe;
mov = addframe(mov,F);
end
%movie(Mf,3)
  mov = close(mov);
  
figure
plot(x,u0,x,u);%axis([0,1,-1.5,1.5])


function A=makematrix(N)
e = ones(N-1,1);
A = spdiags([-e 2*e -e], -1:1, N-1, N-1);

function y=ix(x)
%
% initial data
y=exp(-80*(x-.5).^2);
