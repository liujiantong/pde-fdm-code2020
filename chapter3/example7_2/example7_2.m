% Example 7.2

nx=17; ny=17; lx=1.; ly=1.; hx=lx/(nx-1); hy=ly/(ny-1);
e=zeros(nx,ny); ue=zeros(nx,ny);
tb=cputime;
[u]=helmgoltz_2d_fft('f_e71',1.,0.,lx,ly,nx,ny,'g1_e71','g2_e71','g3_e71','g4_e71');
td=cputime-tb;

% Exact solution
for n=1:nx
    x=(n-1)*hx; xe(n)=x;
    for m=1:ny
        y=(m-1)*hy; ye(m)=y; ue(n,m)=x*(1.-x)*log(1.+y);
        e(n,m)=u(n,m)-ue(n,m);
    end
end
err=norm(e,'fro')/norm(ue,'fro');
fprintf(' Example 7.2  \n');
fprintf('  Relative error - %f   \n',err);
fprintf('  CPU time - %f  \n',td); fprintf(' \n');
mesh(xe,ye,u');
xlabel(' x '); ylabel(' y '); zlabel(' u ');
title(' Example 7.2 ');
clear all;