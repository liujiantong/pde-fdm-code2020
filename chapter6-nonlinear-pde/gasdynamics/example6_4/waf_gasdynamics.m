% Method WAF

function [d,u,p,e,time]=waf_gasdynamics(dd,ud,pd,h,nx,r,cad,tp,gamma,bcl,bcr)

g(1)=0.5*(gamma-1.)/gamma; g(2)=0.5*(gamma+1.)/gamma; g(3)=2.*gamma/(gamma-1.);
g(4)=2./(gamma-1.); g(5)=2./(gamma+1.); g(6)=(gamma-1.)/(gamma+1.);
g(7)=0.5*(gamma-1.); g(8)=1./gamma; g(9)=gamma-1.;
sd(1,1:nx-1)=dd(1:nx-1);
sd(2,1:nx-1)=dd(1:nx-1).*ud(1:nx-1);
sd(3,1:nx-1)=pd(1:nx-1)/(gamma-1.)+0.5*(dd(1:nx-1).*ud(1:nx-1)).*ud(1:nx-1);
av(1:3,1:nx-1)=zeros(3,nx-1);

time=0.;
while time <= tp
    sta=max(abs(sd(2,:)./sd(1,:))+...
        sqrt(gamma*(gamma-1)*(sd(3,:)./sd(1,:)-0.5*(sd(2,:).*sd(2,:))./(sd(1,:).*sd(1,:)))));
    tau=r*sqrt(1.-2.*cad)*h/sta; th=tau/h;
    time=time+tau; %fprintf('  %f  \n',time);
    if bcl == 0
        [f1,f2,f3]=flux_waf(sd(1,1),-sd(2,1),sd(3,1),sd(1,1),sd(2,1),sd(3,1),g(:),th);
    else
        [f1,f2,f3]=flux_waf(sd(1,1),sd(2,1),sd(3,1),sd(1,1),sd(2,1),sd(3,1),g(:),th);
    end    
    fleft(1)=f1; fleft(2)=f2; fleft(3)=f3;
    for n=1:nx-1
        if (n == 1)|(n == nx-1)
            av(:,n)=[ 0.; 0.; 0.;];
        else
            av(:,n)=cad*(sd(:,n+1)-2.*sd(:,n)+sd(:,n-1));
        end    
    end
    for n=1:nx-2
        [f1,f2,f3]=flux_waf(sd(1,n),sd(2,n),sd(3,n),sd(1,n+1),sd(2,n+1),sd(3,n+1),g(:),th);
        fright(1)=f1; fright(2)=f2; fright(3)=f3;
        su(:,n)=sd(:,n)-(fright(:)-fleft(:))*tau/h+av(:,n);
        fleft(:)=fright(:);
    end
    if bcr == 0
        [f1,f2,f3]=flux_waf(sd(1,nx-1),sd(2,nx-1),sd(3,nx-1),sd(1,nx-1),-sd(2,nx-1),...
                                    sd(3,nx-1),g(:),th);
    else
        [f1,f2,f3]=flux_waf(sd(1,nx-1),sd(2,nx-1),sd(3,nx-1),sd(1,nx-1),sd(2,nx-1),...
                                    sd(3,nx-1),g(:),th);
    end    
    fright(1)=f1; fright(2)=f2; fright(3)=f3;    
    su(:,nx-1)=sd(:,nx-1)-(fright(:)-fleft(:))*tau/h;    
    sd(:,1:nx-1)=su(:,1:nx-1);
end
d(1:nx-1)=su(1,1:nx-1); u(1:nx-1)=su(2,1:nx-1)./su(1,1:nx-1);
p(1:nx-1)=(gamma-1.)*(su(3,1:nx-1)-0.5*(su(2,1:nx-1).*su(2,1:nx-1))./su(1,1:nx-1));
e(1:nx-1)=(su(3,1:nx-1)-0.5*(su(2,1:nx-1).*su(2,1:nx-1))./su(1,1:nx-1))./su(1,1:nx-1);
return;