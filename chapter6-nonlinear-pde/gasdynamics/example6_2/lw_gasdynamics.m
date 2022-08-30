% Lax-Wendroff scheme (gasdynamics)

function [d,u,p,e,time]=lw_gasdynamics(dd,ud,pd,h,nx,r,cad,tp,gamma,bcl,bcr)

sd(1,1:nx)=dd(1:nx);
sd(2,1:nx)=dd(1:nx).*ud(1:nx);
sd(3,1:nx)=pd(1:nx)/(gamma-1.)+0.5*(dd(1:nx).*ud(1:nx)).*ud(1:nx);
si(1:3,1:nx+1)=zeros(3,nx+1); av(1:3,1:nx)=zeros(3,nx);
time=0.;
while time <= tp
    sta=max(abs(sd(2,:)./sd(1,:))+...
        sqrt(gamma*(gamma-1)*(sd(3,:)./sd(1,:)-0.5*(sd(2,:).*sd(2,:))./(sd(1,:).*sd(1,:)))));
    tau=r*h*sqrt(1.-2.*cad)/sta;
    time=time+tau; %fprintf('  time - %f  \n',time);
    
% first step
    if bcl == 0
        sl1=sd(1,2); sl2=-sd(2,2); sl3=sd(3,2);
    else
        sl1=sd(1,2); sl2=sd(2,2); sl3=sd(3,2);
    end    
    gl(1)=sl2;
    gl(2)=(gamma-1.)*sl3+0.5*(3.-gamma)*sl2^2/sl1;
    gl(3)=sl2*(gamma*sl3-0.5*(gamma-1.)*sl2^2/sl1)/sl1;
    for n=1:nx
        if (n == 1)|(n == nx)
            av(:,n)=[ 0.; 0.; 0.;];
        else
            av(:,n)=cad*(sd(:,n+1)-2.*sd(:,n)+sd(:,n-1));
        end    
        sr1=sd(1,n); sr2=sd(2,n); sr3=sd(3,n);
        gr(1)=sr2;
        gr(2)=(gamma-1.)*sr3+0.5*(3.-gamma)*sr2^2/sr1;
        gr(3)=sr2*(gamma*sr3-0.5*(gamma-1.)*sr2^2/sr1)/sr1;
        si(1,n)=0.5*(sl1+sr1)-(gr(1)-gl(1))*tau*0.5/h;
        si(2,n)=0.5*(sl2+sr2)-(gr(2)-gl(2))*tau*0.5/h;
        si(3,n)=0.5*(sl3+sr3)-(gr(3)-gl(3))*tau*0.5/h;
        gl(:)=gr(:);
        sl1=sr1; sl2=sr2; sl3=sr3;
    end
    if bcr == 0
        sr1=sd(1,nx-1); sr2=-sd(2,nx-1); sr3=sd(3,nx-1);
    else
        sr1=sd(1,nx-1); sr2=sd(2,nx-1); sr3=sd(3,nx-1);
    end    
    gr(1)=sr2;
    gr(2)=(gamma-1.)*sr3+0.5*(3.-gamma)*sr2^2/sr1;
    gr(3)=sr2*(gamma*sr3-0.5*(gamma-1.)*sr2^2/sr1)/sr1;
    si(1,nx+1)=0.5*(sl1+sr1)-(gr(1)-gl(1))*tau*0.5/h;
    si(2,nx+1)=0.5*(sl2+sr2)-(gr(2)-gl(2))*tau*0.5/h;
    si(3,nx+1)=0.5*(sl3+sr3)-(gr(3)-gl(3))*tau*0.5/h;
    
% second step
    fl(1)=si(2,1);
    fl(2)=(gamma-1.)*si(3,1)+0.5*(3.-gamma)*si(2,1)^2/si(1,1);
    fl(3)=si(2,1)*(gamma*si(3,1)-0.5*(gamma-1.)*si(2,1)^2/si(1,1))/si(1,1);
    for n=2:nx+1
        fr(1)=si(2,n);
        fr(2)=(gamma-1.)*si(3,n)+0.5*(3.-gamma)*si(2,n)^2/si(1,n);
        fr(3)=si(2,n)*(gamma*si(3,n)-0.5*(gamma-1.)*si(2,n)^2/si(1,n))/si(1,n);
        su(:,n-1)=sd(:,n-1)-(fr(:)-fl(:))*tau/h+av(:,n-1);
        fl(:)=fr(:);
    end
    sd(:,1:nx)=su(:,1:nx);
end
d(1:nx)=su(1,1:nx); u(1:nx)=su(2,1:nx)./su(1,1:nx);
p(1:nx)=(gamma-1.)*(su(3,1:nx)-0.5*(su(2,1:nx).*su(2,1:nx))./su(1,1:nx));
e(1:nx)=(su(3,1:nx)-0.5*(su(2,1:nx).*su(2,1:nx))./su(1,1:nx))./su(1,1:nx);
return;