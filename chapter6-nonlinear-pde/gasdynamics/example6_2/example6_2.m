% Example 6.2

xmin=-1.; xmax=1.; nx=101; h=(xmax-xmin)/(nx-1);
r=0.9; cad=0.06; time=0.; gamma=1.4; bcl=1.; bcr=1.;
fprintf(' Example 6.2  \n'); fprintf('  Space step - %f  \n',h); fprintf(' \n');
dd(1:nx)=zeros(1,nx); ud(1:nx)=zeros(1,nx); pd(1:nx)=zeros(1,nx);
x(1:nx)=xmin+((1:nx)-1)*h;

% Test 1
te=' ( Test 1 )';
rol=1.;       ul=0.; pl=1.;
ror=0.125; ur=0.; pr=0.1;
tp=0.5;

% Test 2
%te=' ( Test 2 )';
%rol=1.; ul=0.; pl=1000.;
%ror=1.; ur=0.; pr=1.;
%tp=0.024;

% Test 3
%te=' ( Test 3 )';
%rol=0.1; ul=0.; pl=1.;
%ror=1.; ur=0.; pr=100.;
%tp=0.05;

% Initial conditions
dd(1:(nx-1)/2+1)=rol; ud(1:(nx-1)/2+1)=ul; pd(1:(nx-1)/2+1)=pl;
dd((nx-1)/2+2:nx)=ror;ud((nx-1)/2+2:nx)=ur; pd((nx-1)/2+2:nx)=pr;

[d,u,p,e,time]=lw_gasdynamics(dd,ud,pd,h,nx,r,cad,tp,gamma,bcl,bcr);

% Exact solution
for n=1:nx
    [density, velocity, pressure]=riemann(x(n),time,rol,ul,pl,ror,ur,pr,gamma);
    de(n)=density; ue(n)=velocity; pe(n)=pressure;
    ee(n)=pressure/((gamma-1.)*density);
end
figure(1);
plot(x,de,'b-',x,d,'b--');
xlabel(' x ');
ylabel(' density ');
title([' Example 6.2 ',te]);
figure(2);
plot(x,ue,'b-',x,u,'b--');
xlabel(' x ');
ylabel(' velocity ');
title([' Example 6.2 ',te]);
figure(3);
plot(x,pe,'b-',x,p,'b--');
xlabel(' x ');
ylabel(' pressure ');
title([' Example 6.2 ',te]);
figure(4);
plot(x,ee,'b-',x,e,'b--');
xlabel(' x ');
ylabel(' internal energy ');
title([' Example 6.2 ',te]);
clear all;