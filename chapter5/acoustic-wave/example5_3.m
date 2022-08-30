% Example 5.3

c=1.; ro=1.;  imp=c*ro; l=1.; nx=21; r=0.66666;
h=l/(nx-1); tau=r*h/c; nts=round(1./tau);
uu(1:nx)=zeros(1,nx); um(1:nx)=zeros(1,nx); ud(1:nx)=zeros(1,nx);
pu(1:nx)=zeros(1,nx); pm(1:nx)=zeros(1,nx); pd(1:nx)=zeros(1,nx);
fprintf(' Example 5.3  \n'); fprintf('  \n');
fprintf('  Space step - %f  \n',h);
fprintf('  Time step  - %f  \n',tau);

xa=((1:nx)-1)*h;
for k=1:nts
    time=k*tau;
    [uu,pu]=roe(ud,pd,um,pm,imp,nx,time,r,0,0,'f1_e51','f2_e51');
    ud(1:nx)=um(1:nx); um(1:nx)=uu(1:nx);
    pd(1:nx)=pm(1:nx); pm(1:nx)=pu(1:nx);
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
title(' Example 5.3 ');
clear all;