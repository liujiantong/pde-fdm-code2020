% Solve the discrete Helmgoltz equation using method SOR

function [u,k]=helmgoltz_2d_sor(f,ac,bc,lx,ly,nx,ny,g1,g2,g3,g4)

hx=lx/(nx-1); hy=ly/(ny-1); hx2=1./(hx*hx); hy2=1./(hy*hy);
maxits=3*max(nx,ny); eps=0.01*(min(hx,hy)^2);
x(1:nx)=((1:nx)-1.)*hx; y(1:ny)=((1:ny)-1.)*hy;
v=ones(nx,ny); z=zeros(nx,ny); u=zeros(nx,ny); fa=zeros(nx,ny);
for n=2:nx-1
    for m=2:ny-1
        d(n,m)=feval(ac,x(n),y(m)-0.5*hy)*hy2; c(n,m)=feval(ac,x(n),y(m)+0.5*hy)*hy2;
        b(n,m)=feval(ac,x(n)-0.5*hx,y(m))*hx2; a(n,m)=feval(ac,x(n)+0.5*hx,y(m))*hx2;
        e(n,m)=-(a(n,m)+b(n,m)+c(n,m)+d(n,m))+feval(bc,x(n),y(m));
        fa(n,m)=feval(f,x(n),y(m));
    end
end
for m=1:ny
    u(1,m)=feval(g1,y(m)); u(nx,m)=feval(g3,y(m));
end
for n=1:nx
    u(n,1)=feval(g2,x(n)); u(n,ny)=feval(g4,x(n));
end
err=1.; k=0.; ds=1.; sp=0.;
while err > eps
    q1=0.; q2=0.;
    if ds > eps
        for n=2:nx-1
            for m=2:ny-1
                q2=q2+v(n,m)*v(n,m);
                z(n,m)=-(a(n,m)*v(n+1,m)+b(n,m)*v(n-1,m)+c(n,m)*v(n,m+1)+...
                              d(n,m)*v(n,m-1))/e(n,m);
                q1=q1+z(n,m)*v(n,m);
            end
        end
        sn=q1/q2; ds=abs(sn-sp); sp=sn; v(:,:)=z(:,:)/norm(z,1);
        sj2=sn*sn; omega=2./(1.+sqrt(1.-sj2));
    end    
    nud=0.;
    for n=2:nx-1
        for m=2:ny-1
            residual=a(n,m)*u(n+1,m)+b(n,m)*u(n-1,m)+c(n,m)*u(n,m+1)+...
                         d(n,m)*u(n,m-1)+e(n,m)*u(n,m)-fa(n,m);
            du=omega*residual/e(n,m); nud=nud+du*du;
            u(n,m)=u(n,m)-du;
        end
    end
    err=sqrt(nud)/norm(u,1);
    if k >= maxits
        return;
    else
         k=k+1;
    end
end
return;