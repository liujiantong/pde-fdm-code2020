% Example 5.2

% 0<x<1
%   t
clear all;format long;

c=1.; ro=1.; l=1.; nx=201; r=0.66666; nts=round((nx-1)/r);
h=l/(nx-1); tau=r*h/c; imp=c*ro;
uu(1:nx)=zeros(1,nx); ud(1:nx)=zeros(1,nx);
pu(1:nx)=zeros(1,nx); pd(1:nx)=zeros(1,nx);
fprintf(' Example 5.2  \n'); fprintf('  \n');
fprintf('  Space step - %f  \n',h);
fprintf('  Time step  - %f  \n',tau);

xa=((1:nx)-1)*h;
t = 0.0;
for k=1:nts
    time=k*tau;
    [uu,pu]=lw_acoustics(ud,pd,imp,nx,time,r,0,0,'f1_e51','f2_e51');
    ud(1:nx)=uu(1:nx); pd(1:nx)=pu(1:nx);
    
    if(k==100)
        
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
plot(xe,ue,'b-',xa,uu,'r--');
xlabel(' x ');
ylabel(' u ');
title(' t=time ');
time 

    end
    
     
    
      if(k==200)
        
        % Exact solution
na=4*(nx-1)+1;
for n=1:na
    x=(n-1)*l/(na-1); xe(n)=x; a=time-x/c;
    ue(n)=0.;
    if a >= 0.
       ue(n)=f1_e51(a);
    end   
end

figure(2);
plot(xe,ue,'b-',xa,uu,'r--');
xlabel(' x ');
ylabel(' u ');
title(' t=time ');

time 
      end
    
     if(k==300)
        
        % Exact solution
na=4*(nx-1)+1;
for n=1:na
    x=(n-1)*l/(na-1); xe(n)=x; a=time-x/c;
    ue(n)=0.;
    if a >= 0.
       ue(n)=f1_e51(a);
    end   
end

figure(3);
plot(xe,ue,'b-',xa,uu,'r--');
xlabel(' x ');
ylabel(' u ');
title(' t=time ');
time 
    end
    
end


