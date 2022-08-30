% Example 5.1

c=1.; ro=1.; l=1.; nx=21; r=0.66666; nts=round((nx-1)/r);
h=l/(nx-1); tau=r*h/c;
ud(1:nx)=zeros(1,nx); um(1:nx)=zeros(1,nx);
uu(1:nx)=zeros(1,nx);
fprintf(' Example 5.1  \n'); fprintf('  \n');
fprintf('  Space step - %f  \n',h);
fprintf('  Time step  - %f  \n',tau);

xa=((1:nx)-1)*h;  
for k=1:nts
    time=k*tau;
    [uu]=cross(ud,um,c,ro,nx,h,time,tau,0,0,'f1_e51','f2_e51');
    ud(1:nx)=um(1:nx); um(1:nx)=uu(1:nx);    
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
title(' Example 5.1 ');
clear all;