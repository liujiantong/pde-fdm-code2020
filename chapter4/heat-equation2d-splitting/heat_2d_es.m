% Two-dimensional equation of heat transfer for uniform medium
% Explicit splitting scheme

function [uu]=heat_2d_es(ud,kappa,a,nx,ny,hx,hy,tau,time,bxl,bxr,fbxl,fbxr,byl,byr,fbyl,fbyr,fs)

x=((1:nx)-1)*hx; y=((1:ny)-1)*hy; v(1:nx,1:ny)=zeros(nx,ny);
bx=kappa*tau/(hx*hx); bxa=2.*bx; by=kappa*tau/(hy*hy); bya=2.*by;
if max(bx,by) > 0.5
    disp('Scheme is unstable'); return;
end
ms=1;
if byl == 0
    ms=2;
end    
me=ny;
if byr == 0
    me=ny-1;
end    
for m=ms:me
    if bxl == 0
        v(1,m)=feval(fbxl,y(m),time);
    else
        v(1,m)=(1.-bxa)*ud(1,m)+bxa*ud(2,m)-bxa*hx*feval(fbxl,y(m),time-tau)/a;
    end    
    uu(1,m)=v(1,m);
    for n=2:nx-1            
        v(n,m)=ud(n,m)+bx*(ud(n+1,m)-2.*ud(n,m)+ud(n-1,m))+tau*feval(fs,x(n),y(m),time-tau);
    end    
    if bxr == 0
        v(nx,m)=feval(fbxr,y(m),time);
    else
        v(nx,m)=(1.-bxa)*ud(nx,m)+bxa*ud(nx-1,m)+bxa*hx*feval(fbxr,y(m),time-tau)/a;
    end
    uu(nx,m)=v(nx,m);
end
ns=1;
if bxl == 0
    ns=2;
end    
ne=nx;
if bxr == 0
    ne=nx-1;
end    
for n=ns:ne    
    if byl == 0
        uu(n,1)=feval(fbyl,x(n),time);
    else
        uu(n,1)=(1.-bya)*v(n,1)+bya*v(n,2)-bya*hy*feval(fbyl,x(n),time-tau)/a;
    end
    for m=2:ny-1        
        uu(n,m)=v(n,m)+by*(v(n,m+1)-2.*v(n,m)+v(n,m-1));
    end
    if byr == 0
        uu(n,ny)=feval(fbyr,x(n),time);
    else
        uu(n,ny)=(1.-bya)*v(n,ny)+bya*v(n,ny-1)+bya*hy*feval(fbyr,x(n),time-tau)/a;
    end
end