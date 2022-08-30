% Solve the Helmgoltz equation using FFT

function [u]=helmgoltz_2d_fft(f,a,b,lx,ly,nx,ny,g1,g2,g3,g4)

hx=lx/(nx-1); hy=ly/(ny-1); hxs=1./(hx*hx); hys=1./(hy*hy);
c1=a*hxs; c2=a*hys; c3=2.*(c1+c2)-b;
u=zeros(nx,ny);
x(1:nx)=((1:nx)-1.)*hx; y(1:ny)=((1:ny)-1.)*hy;
d=zeros(nx-2,ny);

for n=1:nx-2
    for m=1:ny
        d(n,m)=feval(f,x(n+1),y(m));
    end
end
for m=1:ny
    d(1,m)=d(1,m)-c2*feval(g1,y(m)); d(nx-2,m)=d(nx-2,m)-c2*feval(g3,y(m));
end
for n=1:nx-2
    d(n,1)=feval(g2,x(n+1));
end
for n=1:nx-2
    d(n,ny)=feval(g4,x(n+1));
end
d = fast_sin_transform(d);
for k=1:nx-2
    p=c1; r=c1; s=-(c3-2.*c2*cos(pi*k/(nx-1)));
    for m=1:ny
        fs(m)=d(k,m);
    end
    fs=sweep(p,s,r,1.,0.,0.,1.,fs);
    d(k,1:ny)=fs(1:ny);
end
d = fast_sin_transform(d); u(2:nx-1,1:ny)=2.*d(1:nx-2,1:ny)/(nx-1);
for m=1:ny
    u(1,m)=feval(g1,y(m)); u(nx,m)=feval(g3,y(m));
end
return;