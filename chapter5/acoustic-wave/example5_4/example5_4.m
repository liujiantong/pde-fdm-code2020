% Example 5.4

c=1.; ro=1.; l=1.; nx=601; r=0.66666; nts=round((nx-1)/r);
h=l/(nx-1); tau=r*h/c; imp=c*ro;
ud(1:nx-1)=zeros(1,nx-1); uu(1:nx-1)=zeros(1,nx-1);
pd(1:nx-1)=zeros(1,nx-1); pu(1:nx-1)=zeros(1,nx-1);
fprintf(' Example 5.4  \n'); fprintf('  \n');
fprintf('  Space step - %f  \n',h);
fprintf('  Time step  - %f  \n',tau);

xa=((1:nx-1)-0.5)*h;
for k=1:nts
    time=k*tau;
    [uu,pu]=godunov_acoustics(ud,pd,imp,nx,time,r,0,0,'f1_e51','f2_e51');   
    ud(1:nx-1)=uu(1:nx-1); pd(1:nx-1)=pu(1:nx-1);
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
title(' Example 5.4 ');
clear all;