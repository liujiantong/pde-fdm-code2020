% Example 5.5

c=1.; ro=1.; l=1.; nx=21; r=0.66666; nts=round((nx-1)/r);
h=l/(nx-1); tau=r*h/c; imp=c*ro;
ua(1:nx-1)=zeros(1,nx-1); ub(1:nx-1)=zeros(1,nx-1);
ud(1:nx)=zeros(1,nx); uu(1:nx)=zeros(1,nx);
pa(1:nx-1)=zeros(1,nx-1); pb(1:nx-1)=zeros(1,nx-1);
pd(1:nx)=zeros(1,nx); pu(1:nx)=zeros(1,nx);
z1=exp(-1.); z2=1.-z1; z3=1.-exp(-r); z4=z1*(exp(r)-1.);
pc=z2*(3.*z1-1.);
a1=(-z2*z2+z1*(1.-z3)+(1.-2.*z1)*(z1+z4))/pc;
a2=z2*(z3-z4)/pc;
a3=(-z2*z2+z1*(z1+z4)+(1.-2.*z1)*(1.-z3))/pc;
b1=h*(z1*z3+(1.-2.*z1)*z4-z2*z2*r)/(c*pc);
b2=h*(r*z2*(1.+z1)-z2*(z3+z4))/(c*pc);
b3=h*(z1*z4+(1.-2.*z1)*z3-z2*z2*r)/(c*pc);
fprintf(' Example 5.5  \n'); fprintf('  \n');
fprintf('  Space step - %f  \n',h);
fprintf('  Time step  - %f  \n',tau);

xa=((1:nx)-1)*h;
for k=1:nts
    time=k*tau;
    [ub,pb,uu,pu]=pem(ua,pa,ud,pd,imp,nx,r,time,tau,a1,a2,a3,b1,b2,b3,...
                                 0,0,'f1_e51','f2_e51');
    ua(1:nx-1)=ub(1:nx-1); pa(1:nx-1)=pb(1:nx-1);
    ud(1:nx)=uu(1:nx); pd(1:nx)=pu(1:nx);
end

% Exact solution
na=4*(nx-1)+1;
for n=1:na
    x=(n-1)*l/(na-1); xe(n)=x; a=time-x/c;
    ue(n)=0.;
    if a >= 0.
       ue(n)=f1_e51(a);
    end   
end

figure(1);
plot(xe,ue,'k-',xa,uu,'k--');
xlabel(' x ');
ylabel(' u ');
title(' Example 5.5 ');
clear all;