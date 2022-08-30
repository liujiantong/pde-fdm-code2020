% Method WAF (TVD)

function [d,u,p,e,time]=waf_tvd_gasdynamics(dd,ud,pd,h,nx,r,tp,gamma,bcl,bcr,limiter)

g(1)=0.5*(gamma-1.)/gamma; g(2)=0.5*(gamma+1.)/gamma; g(3)=2.*gamma/(gamma-1.);
g(4)=2./(gamma-1.); g(5)=2./(gamma+1.); g(6)=(gamma-1.)/(gamma+1.);
g(7)=0.5*(gamma-1.); g(8)=1./gamma; g(9)=gamma-1.;
sd(1,1:nx-1)=dd(1:nx-1);
sd(2,1:nx-1)=dd(1:nx-1).*ud(1:nx-1);
sd(3,1:nx-1)=pd(1:nx-1)/(gamma-1.)+0.5*(dd(1:nx-1).*ud(1:nx-1)).*ud(1:nx-1);
rn(1:4,1:nx-2)=zeros(4,nx-2); un(1:4,1:nx-2)=zeros(4,nx-2); pn(1:4,1:nx-2)=zeros(4,nx-2);
cp(1:3,1:nx-2)=zeros(3,nx-2);

time=0.;
while time <= tp
    sta=max(abs(sd(2,:)./sd(1,:))+...
        sqrt(gamma*(gamma-1)*(sd(3,:)./sd(1,:)-0.5*(sd(2,:).*sd(2,:))./(sd(1,:).*sd(1,:)))));
    tau=r*h/sta; th=tau/h;
    time=time+tau; %fprintf('  %f  \n',time);
    for n=1:nx-2
        [rw,uw,pw,cw]=flux_waf_tvd(sd(1,n),sd(2,n),sd(3,n),sd(1,n+1),sd(2,n+1),sd(3,n+1),g(:));
        rn(:,n)=rw(:); un(:,n)=uw(:); pn(:,n)=pw(:); cp(:,n)=cw(:);
    end    
    if bcl == 0
        [f1,f2,f3]=flux_godunov(sd(1,1),-sd(2,1),sd(3,1),sd(1,1),sd(2,1),sd(3,1),g(:));
    else
        [f1,f2,f3]=flux_godunov(sd(1,1),sd(2,1),sd(3,1),sd(1,1),sd(2,1),sd(3,1),g(:));
    end    
    fleft(1)=f1; fleft(2)=f2; fleft(3)=f3;
    for n=1:nx-2
        sum1=rn(1,n)+rn(4,n); sum2=un(1,n)+un(4,n); sum3=pn(1,n)+pn(4,n);
        for p=1:3
            ddl=rn(p+1,n)-rn(p,n);
            if abs(ddl) < 1.0e-6
                if ddl < 0.
                    ddl=-1.0e-6;
                else
                    ddl=1.0e-6;
                end    
            end    
            udl=un(p+1,n)-un(p,n);
            if abs(udl) < 1.0e-6
                if udl < 0.
                    udl=-1.0e-6;
                else
                    udl=1.0e-6;
                end    
            end    
            pdl=pn(p+1,n)-pn(p,n);
            if abs(pdl) < 1.0e-6
                if pdl < 0.
                    pdl=-1.0e-6;
                else
                    pdl=1.0e-6;
                end    
            end    
            if cp(p,n) >= 0.
                if n==1
                    q1=0.; q2=0.; q3=0.;
                else
                    q1=(rn(p+1,n-1)-rn(p,n-1))/ddl; q2=(un(p+1,n-1)-un(p,n-1))/udl;
                    q3=(pn(p+1,n-1)-pn(p,n-1))/pdl;
                end
                b=cp(p,n)*tau/h;
                fi1=feval(limiter,q1,b); fi2=feval(limiter,q2,b); fi3=feval(limiter,q3,b);
                sum1=sum1+fi1*(rn(p,n)-rn(p+1,n)); sum2=sum2+fi2*(un(p,n)-un(p+1,n));
                sum3=sum3+fi3*(pn(p,n)-pn(p+1,n));
            else
                if n==nx-2
                    q1=0.; q2=0.; q3=0.;
                else
                    q1=(rn(p+1,n+1)-rn(p,n+1))/ddl; q2=(un(p+1,n+1)-un(p,n+1))/udl;
                    q3=(pn(p+1,n+1)-pn(p,n+1))/pdl;
                end
                b=abs(cp(p,n))*tau/h;
                fi1=feval(limiter,q1,b); fi2=feval(limiter,q2,b); fi3=feval(limiter,q3,b);
                sum1=sum1-fi1*(rn(p,n)-rn(p+1,n)); sum2=sum2-fi2*(un(p,n)-un(p+1,n));
                sum3=sum3-fi3*(pn(p,n)-pn(p+1,n));
            end    
        end
        rf=0.5*sum1; uf=0.5*sum2; pf=0.5*sum3;
        fright(1)=rf*uf; fright(2)=pf+rf*uf*uf; fright(3)=0.5*uf*(pf*g(3)+rf*uf*uf);
        su(:,n)=sd(:,n)-(fright(:)-fleft(:))*tau/h;
        fleft(:)=fright(:);
    end
    if bcr == 0
        [f1,f2,f3]=flux_godunov(sd(1,nx-1),sd(2,nx-1),sd(3,nx-1),sd(1,nx-1),-sd(2,nx-1),sd(3,nx-1),g(:));
    else
        [f1,f2,f3]=flux_godunov(sd(1,nx-1),sd(2,nx-1),sd(3,nx-1),sd(1,nx-1),sd(2,nx-1),sd(3,nx-1),g(:));
    end    
    fright(1)=f1; fright(2)=f2; fright(3)=f3;    
    su(:,nx-1)=sd(:,nx-1)-(fright(:)-fleft(:))*tau/h;    
    sd(:,1:nx-1)=su(:,1:nx-1);
end
d(1:nx-1)=su(1,1:nx-1); u(1:nx-1)=su(2,1:nx-1)./su(1,1:nx-1);
p(1:nx-1)=(gamma-1.)*(su(3,1:nx-1)-0.5*(su(2,1:nx-1).*su(2,1:nx-1))./su(1,1:nx-1));
e(1:nx-1)=(su(3,1:nx-1)-0.5*(su(2,1:nx-1).*su(2,1:nx-1))./su(1,1:nx-1))./su(1,1:nx-1);
return;