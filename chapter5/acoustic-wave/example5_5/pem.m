% PEM scheme

function [ub,pb,uu,pu]=pem(ua,pa,ud,pd,imp,nx,r,time,tau,a1,a2,a3,b1,b2,b3,bl,br,fbl,fbr)

if r > 1.
    disp('Scheme is unstable'); fprintf('  \n'); return;
end
rpa=pd(1)+imp*ud(1); rpb=pa(1)+imp*ua(1); rpc=pd(2)+imp*ud(2);
if bl == 0
    uu(1)=feval(fbl,time); 
    ful=(feval(fbl,time)+feval(fbl,time-tau))*0.5*tau;
    pu(1)=a1*rpa+a2*rpb+a3*rpc-imp*uu(1);
    fpl=b1*rpa+b2*rpb+b3*rpc-imp*ful;
else
    pu(1)=feval(fbl,time); 
    fpl=(feval(fbl,time)+feval(fbl,time-tau))*0.5*tau;
    uu(1)=(a1*rpa+a2*rpb+a3*rpc-pu(1))/imp;
    ful=(b1*rpa+b2*rpb+b3*rpc-fpl)/imp;
end
for n=1:nx-1
    if n == nx-1
        rma=pd(n)-imp*ud(n); rmb=pa(n)-imp*ua(n); rmc=pd(n+1)-imp*ud(n+1);
        if br == 0
            uu(nx)=feval(fbr,time);
            fur=(feval(fbr,time)+feval(fbr,time-tau))*0.5*tau;
            pu(nx)=a3*rma+a2*rmb+a1*rmc+imp*uu(nx);
            fpr=b3*rma+b2*rmb+b1*rmc+imp*fur;
        else
            pu(nx)=feval(fbr,time);
            fpr=(feval(fbr,time)+feval(fbr,time-tau))*0.5*tau;
            uu(nx)=(pu(nx)-a3*rma-a2*rmb-a1*rmc)/imp;
            fur=(fpr-b3*rma-b2*rmb-b1*rmc)/imp;
        end
    else
        rpa=pd(n+1)+imp*ud(n+1); rpb=pa(n+1)+imp*ua(n+1); rpc=pd(n+2)+imp*ud(n+2);
	    rma=pd(n)-imp*ud(n); rmb=pa(n)-imp*ua(n); rmc=pd(n+1)-imp*ud(n+1);
	    pu(n+1)=(a1*rpa+a2*rpb+a3*rpc+a3*rma+a2*rmb+a1*rmc)*0.5;
	    uu(n+1)=(a1*rpa+a2*rpb+a3*rpc-a3*rma-a2*rmb-a1*rmc)*0.5/imp;
	    fpr=(b1*rpa+b2*rpb+b3*rpc+b3*rma+b2*rmb+b1*rmc)*0.5;
        fur=(b1*rpa+b2*rpb+b3*rpc-b3*rma-b2*rmb-b1*rmc)*0.5/imp;
    end
    ub(n)=ua(n)+(fpr-fpl)*r/(imp*tau); pb(n)=pa(n)+(fur-ful)*r*imp/tau;
    ful=fur; fpl=fpr;
end
return;