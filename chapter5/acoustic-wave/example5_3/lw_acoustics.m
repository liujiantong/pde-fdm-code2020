% Lax-Wendroff scheme (acoustics)

function [uu,pu]=lw_acoustics(ud,pd,imp,nx,time,r,bl,br,fbl,fbr)

if r > 1
    disp('Scheme is unstable'); fprintf('  \n'); return;
end
if bl == 0       
    uu(1)=feval(fbl,time);
    pu(1)=pd(1)-(ud(2)-feval(fbl,time))*imp*r+(pd(2)-pd(1))*r*r;
else
    pu(1)=feval(fbl,time);
    uu(1)=ud(1)-(pd(2)-feval(fbl,time))*r/imp+(ud(2)-ud(1))*r*r;
end
for n=2:nx-1
    pu(n)=pd(n)-(ud(n+1)-ud(n-1))*imp*r*0.5+(pd(n+1)-2.*pd(n)+pd(n-1))*r*r*0.5;
    uu(n)=ud(n)-(pd(n+1)-pd(n-1))*r*0.5/imp+ (ud(n+1)-2.*ud(n)+ud(n-1))*r*r*0.5;
end
if br == 0
    uu(nx)=feval(fbr,time);
    pu(nx)=pd(nx)-(feval(fbr,time)-ud(nx-1))*imp*r+(pd(nx-1)-pd(nx))*r*r;
else
    pu(nx)=feval(fbr,time);
    uu(nx)=ud(nx)-(feval(fbr,time)-pd(nx-1))*r/imp+(ud(nx-1)-ud(nx))*r*r;
end
return;