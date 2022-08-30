% Solve the Poisson equation using FFT (for ns_2d_vsf.m)

function [u]=poisson_ns_fft(f,nx,ny,hx,hy)

c1=1./(hx*hx); c2=1./(hy*hy); c3=2.*(c1+c2);
u=zeros(nx,ny); d=zeros(nx-2,ny);
for n=1:nx-2
    for m=1:ny
        d(n,m)=-f(n+1,m);
    end
end
for n=1:nx-2
    d(n,1)=0.;
end
for n=1:nx-2
    d(n,ny)=0.;
end
d = fast_sin_transform(d);
for k=1:nx-2
    p=c1; r=c1; s=-(c3-2.*c2*cos(pi*k/(nx-1)));
    for m=1:ny
        fs(m)=d(k,m);
    end
    gamma1=1.; gamma2=0.; gamma4=1.; gamma3=0.;
    fs=sweep(p,s,r,gamma1,gamma2,gamma3,gamma4,fs);
    d(k,1:ny)=fs(1:ny);
end
d = fast_sin_transform(d); u(2:nx-1,1:ny)=2.*d(1:nx-2,1:ny)/(nx-1);
for m=1:ny
    u(1,m)=0.; u(nx,m)=0.;
end
return;