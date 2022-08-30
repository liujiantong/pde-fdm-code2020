function wave2d(N,M,T)
%
% Solve the wave equation on a periodic 2D domain
% using leapfrog
dt=1/M;
h=1/N;
disp(sprintf('CFL parameter dt/h = %g',dt/h))
x=linspace(0,1,N+1)';
[X,Y]=meshgrid(x,x);
u0=ic(X,Y);
u=u0;
uprev=u0;
unew=u0;
%mesh(X,Y,u0);
mesh(X,Y,u0);axis('square'),colorbar;axis([0,1,0,1,-1,1]);drawnow;pause
%Mf(1)=getframe;
v=zeros(N,1);
for j=1:ceil(M*T)
for r=2:N
unew(r,2:N)=2*u(r,2:N)-uprev(r,2:N)+...
dt^2*(u(r+1,2:N)+u(r,(2:N)+1)+u(r,(2:N)-1)+u(r-1,2:N)-4*u(r,2:N))/h^2;
end
uprev=u;
u=unew;
mesh(X,Y,u);axis('square');colorbar,axis([0,1,0,1,-1,1]);drawnow;%pause
end
%movie(Mf,3)
figure
mesh(X,Y,u);axis('square');axis([0,1,0,1,-1,1]);colorbar
figure
contourf(X,Y,u);axis('square');colorbar
norm(u-u0,'fro')/norm(u0,'fro')
function y=ic(X,Y)
%y=exp(-80*((X-.5).^2+(Y-0.5).^2));
y=sin(2*pi*X).*sin(2*pi*Y);
