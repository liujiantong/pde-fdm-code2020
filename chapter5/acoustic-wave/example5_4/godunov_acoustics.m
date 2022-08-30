% Scheme of S.K.Godunov of the 1st order (acoustics)

function [uu,pu]=godunov_acoustics(ud,pd,imp,nx,time,r,bl,br,fbl,fbr)

if r > 1.
    disp('Scheme is unstable'); fprintf('  \n'); return;
end    
rm=pd(1)-imp*ud(1);
if bl == 0
    ul=feval(fbl,time); pl=rm+imp*ul;
else
    pl=feval(fbl,time); ul=(pl-rm)/imp;
end
for n=1:nx-2
    ur=0.5*(ud(n)+ud(n+1))+0.5*(pd(n)-pd(n+1))/imp;
    pr=0.5*(ud(n)-ud(n+1))*imp+0.5*(pd(n)+pd(n+1))/imp;
    uu(n)=ud(n)-(pr-pl)*r/imp; pu(n)=pd(n)-(ur-ul)*r*imp;
    ul=ur; pl=pr;
end
rp=pd(nx-1)+imp*ud(nx-1);
if br == 0
    ur=feval(fbr,time); pr=rp-imp*ur;
else
    pr=feval(fbr,time); ur=(rp-pr)/imp;
end
uu(nx-1)=ud(nx-1)-(pr-pl)*r/imp;
pu(nx-1)=pd(nx-1)-(ur-ul)*r*imp;
return;