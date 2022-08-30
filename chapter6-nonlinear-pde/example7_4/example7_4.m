% Example 7.4

nx=11; ny=11; lx=1.; ly=1.; hx=lx/(nx-1); hy=ly/(ny-1);
ui=zeros(nx,ny); ue=zeros(nx,ny); e=zeros(nx,ny);

[u,k]=poisson_2d_td(ui,'f_e74','a_e74',lx,ly,nx,ny,'g1_e74','g2_e74','g3_e74','g4_e74');

% Exact solution
for n=1:nx
    x=(n-1)*hx; xe(n)=x;
    for m=1:ny
        y=(m-1)*hy; ye(m)=y; ue(n,m)=log(1.+x*y);
        e(n,m)=u(n,m)-ue(n,m);
    end
end
err=norm(e,'fro')/norm(ue,'fro');
fprintf(' Example 7.4  \n');
fprintf('  Relative error - %f   \n',err);
fprintf('  Number of iterations - %i  \n',k); fprintf(' \n');
figure;mesh(xe,ye,u');

xlabel(' x '); ylabel(' y '); zlabel(' u ');
title(' Example 7.4 ');

figure;mesh(xe,ye,ue');
xlabel(' x '); ylabel(' y '); zlabel(' uexact ');

clear all;