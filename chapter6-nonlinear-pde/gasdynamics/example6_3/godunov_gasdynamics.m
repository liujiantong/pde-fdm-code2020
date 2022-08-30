% Method of S. K. Godunov

function [d,u,p,e,time]=godunov_gasdynamics(dd,ud,pd,h,nx,r,tp,gamma,bcl,bcr)

g(1)=0.5*(gamma-1.)/gamma; g(2)=0.5*(gamma+1.)/gamma; g(3)=2.*gamma/(gamma-1.);
g(4)=2./(gamma-1.); g(5)=2./(gamma+1.); g(6)=(gamma-1.)/(gamma+1.);
g(7)=0.5*(gamma-1.); g(8)=1./gamma; g(9)=gamma-1.;
sd(1,1:nx-1)=dd(1:nx-1);
sd(2,1:nx-1)=dd(1:nx-1).*ud(1:nx-1);
sd(3,1:nx-1)=pd(1:nx-1)/(gamma-1.)+0.5*(dd(1:nx-1).*ud(1:nx-1)).*ud(1:nx-1);

time=0.;
while time <= tp
    sta=max(abs(sd(2,:)./sd(1,:))+...
        sqrt(gamma*(gamma-1)*(sd(3,:)./sd(1,:)-0.5*(sd(2,:).*sd(2,:))./(sd(1,:).*sd(1,:)))));
    tau=r*h/sta;
    time=time+tau; %fprintf('  %f  \n',time);
    if bcl == 0
        [f1,f2,f3]=flux_godunov(sd(1,1),-sd(2,1),sd(3,1),sd(1,1),sd(2,1),sd(3,1),g(:));
    else
        [f1,f2,f3]=flux_godunov(sd(1,1),sd(2,1),sd(3,1),sd(1,1),sd(2,1),sd(3,1),g(:));
    end    
    fleft(1)=f1; fleft(2)=f2; fleft(3)=f3;
    for n=1:nx-2
        [f1,f2,f3]=flux_godunov(sd(1,n),sd(2,n),sd(3,n),sd(1,n+1),sd(2,n+1),sd(3,n+1),g(:));
        fright(1)=f1; fright(2)=f2; fright(3)=f3;
        su(:,n)=sd(:,n)-(fright(:)-fleft(:))*tau/h;
        fleft(:)=fright(:);
    end
    if bcr == 0
        [f1,f2,f3]=flux_godunov(sd(1,nx-1),sd(2,nx-1),sd(3,nx-1),sd(1,nx-1),-sd(2,nx-1),...
                                            sd(3,nx-1),g(:));
    else
        [f1,f2,f3]=flux_godunov(sd(1,nx-1),sd(2,nx-1),sd(3,nx-1),sd(1,nx-1),sd(2,nx-1),...
                                            sd(3,nx-1),g(:));
    end    
    fright(1)=f1; fright(2)=f2; fright(3)=f3;    
    su(:,nx-1)=sd(:,nx-1)-(fright(:)-fleft(:))*tau/h;    
    sd(:,1:nx-1)=su(:,1:nx-1);
end
d(1:nx-1)=su(1,1:nx-1); u(1:nx-1)=su(2,1:nx-1)./su(1,1:nx-1);
p(1:nx-1)=(gamma-1.)*(su(3,1:nx-1)-0.5*(su(2,1:nx-1).*su(2,1:nx-1))./su(1,1:nx-1));
e(1:nx-1)=(su(3,1:nx-1)-0.5*(su(2,1:nx-1).*su(2,1:nx-1))./su(1,1:nx-1))./su(1,1:nx-1);
return;