% Sweep method for tridiagonal matrix

function u=sweep(p,s,r,gamma1,gamma2,gamma3,gamma4,d)

nx= length(d);
alfa(1:nx)=zeros(1,nx); beta(1:nx)=zeros(1,nx);
alfa(1)=-gamma2/gamma1; beta(1)=d(1)/gamma1;
for n=2:nx-1
    aa=1./(s+p*alfa(n-1));
    alfa(n)=-r*aa; beta(n)=(d(n)-p*beta(n-1))*aa;
end
aa=1./(gamma4+gamma3*alfa(nx-1));
beta(nx)=(d(nx)-gamma3*beta(nx-1))*aa;
u(nx)=beta(nx);
for n=nx-1:-1:1
    u(n)=alfa(n)*u(n+1)+beta(n);
end
return;