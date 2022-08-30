% Solve the nonlinear Poisson equation using the method of time development

function [u,k]=poisson_2d_td(up,f,a,lx,ly,nx,ny,g1,g2,g3,g4)

hx=lx/(nx-1); hy=ly/(ny-1); hxs=1./(hx*hx); hys=1./(hy*hy);
dx(1:nx)=zeros(1,nx); dy(1:ny)=zeros(1,ny);
px(1:nx)=zeros(1,nx); rx(1:nx)=zeros(1,nx); sx(1:nx)=zeros(1,nx);
py(1:ny)=zeros(1,ny); ry(1:ny)=zeros(1,ny); sy(1:ny)=zeros(1,ny);
x(1:nx)=((1:nx)-1.)*hx; y(1:ny)=((1:ny)-1.)*hy;

u(:,:)=up(:,:);
for n=1:nx
    for m=1:ny
        fs(n,m)=feval(f,x(n),y(m));
    end
end
for m=1:ny
    u(1,m)=feval(g1,y(m)); u(nx,m)=feval(g3,y(m));
end
for n=1:nx
    u(n,1)=feval(g2,x(n)); u(n,ny)=feval(g4,x(n));
end
dif=1.; k=0;
while dif > 1.e-5
    k=k+1;
    for n=1:nx
        for m=1:ny
            as(n,m)=feval(a,u(n,m),x(n),y(m));
        end
    end
    tau=(nx+ny)*hx*hy/(pi*max(max(as))+min(min(as)));
    v(1:nx,1)=u(1:nx,1);
    for m=2:ny-1
        dx(1)=u(1,m);
        for n=2:nx-1
            px(n)=-0.25*tau*(as(n,m)+as(n-1,m))*hxs;
            rx(n)=-0.25*tau*(as(n,m)+as(n+1,m))*hxs;
            sx(n)=1.-px(n)-rx(n);
            dx(n)=u(n,m)+0.25*tau*(as(n,m)+as(n,m+1))*(u(n,m+1)-u(n,m))*hys-...
                                  0.25*tau*(as(n,m)+as(n,m-1))*(u(n,m)-u(n,m-1))*hys-0.5*tau*fs(n,m);
        end
        dx(nx)=u(nx,m);
        dx=sweepa(px,sx,rx,1.,0.,0.,1.,dx);
        v(1:nx,m)=dx';
    end
    v(1:nx,ny)=u(1:nx,ny);
    for n=2:nx-1
        dy(1)=u(n,1);
        for m=2:ny-1
            py(m)=-0.25*tau*(as(n,m)+as(n,m-1))*hys;
            ry(m)=-0.25*tau*(as(n,m)+as(n,m+1))*hys;
            sy(m)=1.-py(m)-ry(m);
            dy(m)=v(n,m)+0.25*tau*(as(n,m)+as(n+1,m))*(v(n+1,m)-v(n,m))*hxs-...
                                  0.25*tau*(as(n,m)+as(n-1,m))*(v(n,m)-v(n-1,m))*hxs-0.5*tau*fs(n,m);
        end
        dy(ny)=u(n,ny);
        dy=sweepa(py,sy,ry,1.,0.,0.,1.,dy);
        u(n,1:ny)=dy(1:ny);
    end
    dif=norm(u-up,1)/norm(u,1); up(:,:)=u(:,:);
end
return;