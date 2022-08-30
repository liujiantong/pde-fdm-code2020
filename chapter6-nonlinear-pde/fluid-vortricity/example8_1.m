% Example 8.1

lx=1.; ly=1.; nx=33; ny=33; hx=lx/(nx-1); hy=ly/(ny-1);
Re=1.0e3; tp=15.;
u(1:nx,1:ny)=zeros(nx,ny); v(1:nx,1:ny)=zeros(nx,ny); w(1:nx,1:ny)=zeros(nx,ny); 
psi(1:nx,1:ny)=zeros(nx,ny);
x(1:nx)=((1:nx)-1.)*hx; y(1:ny)=((1:ny)-1.)*hy;
fprintf(' Example 8.1  \n');
fprintf('  Space step in the x direction - %f  \n',hx);
fprintf('  Space step in the y direction - %f  \n',hy);  fprintf('   \n');
fprintf('  Please wait...computations are in progress \n');
fprintf('  ..... \n');

[u,v,w,psi,time]=ns_2d_vsf(u,v,w,Re,tp,lx,ly,nx,ny,1,1,1,0,'g1_e81','g2_e81','g3_e81','g4_e81');

%figure(1);
%contour(x,y,psi',10,'-k');
xlabel(' x ');
ylabel(' y ');
title(' Example 8.1 ');
clear all;