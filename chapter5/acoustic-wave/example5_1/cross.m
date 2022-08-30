% Scheme 'cross'

function [uu]=cross(ud,um,c,ro,nx,h,time,tau,bl,br,fbl,fbr)

r=c*tau/h;
if r > 1
    disp('Scheme is unstable'); fprintf('  \n'); return;
end
if bl == 0
    um(1)=feval(fbl,time-tau); uu(1)=feval(fbl,time);
else
    uu(1)=2.*um(1)-ud(1)+2.*r*r*(um(2)-um(1)-h*feval(fbl,time-tau)/(c*c*ro));
end
for n=2:nx-1
    uu(n)=2.*um(n)-ud(n)+r*r*(um(n+1)-2.*um(n)+um(n-1));
end
if br == 0
    um(nx)=feval(fbr,time-tau); uu(nx)=feval(fbr,time);
else
    uu(nx)=2.*um(nx)-ud(nx)+...
                2.*r*r*(um(nx-1)-um(nx)+h*feval(fbr,time-tau)/(c*c*ro));
end