% Scheme of P.Roe

function [uu,pu]=roe(ud,pd,um,pm,imp,nx,time,r,bl,br,fbl,fbr)

if r > 1.
    disp('Scheme is unstable'); fprintf('  \n'); return;
end
ra=1.-2.*r;
if bl == 0
    uu(1)=feval(fbl,time);
    rma=pm(1)-imp*um(1); rmb=pm(2)-imp*um(2); rmc=pd(2)-imp*ud(2);
    pu(1)=imp*uu(1)+rmc+(rma-rmb)*ra;
else
    pu(1)=feval(fbl,time);
    rma=pm(1)-imp*um(1); rmb=pm(2)-imp*um(2); rmc=pd(2)-imp*ud(2);
    uu(1)=(pu(1)-rmc-(rma-rmb)*ra)/imp;
end
for n=2:nx-1
    rpa=pm(n)+imp*um(n); rpb=pm(n-1)+imp*um(n-1); rpc=pd(n-1)+imp*ud(n-1);
    rma=pm(n)-imp*um(n); rmb=pm(n+1)-imp*um(n+1); rmc=pd(n+1)-imp*ud(n+1);
    pu(n)=(rpc+rmc)*0.5+(rpa+rma-rpb-rmb)*ra*0.5;
    uu(n)=(rpc-rmc)*0.5/imp+(rpa-rma-rpb+rmb)*ra*0.5/imp;
end
if br == 0
    uu(nx)=feval(fbr,time);
    rpa=pm(nx)+imp*um(nx); rpb=pm(nx-1)+imp*um(nx-1); rpc=pd(nx-1)+imp*ud(nx-1);
    pu(nx)=rmc+(rpa-rpb)*ra-imp*uu(nx);
else
    pu(nx)=feval(fbr,time);
    rpa=pm(nx)+imp*um(nx); rpb=pm(nx-1)+imp*um(nx-1); rpc=pd(nx-1)+imp*ud(nx-1);
    uu(nx)=(-pu(nx)+rpc+(rpa-rpb)*ra)/imp;
end
return;