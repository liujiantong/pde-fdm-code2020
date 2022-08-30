% Example 6.5

xmin=-1.; xmax=1.; nx=101; h=(xmax-xmin)/(nx-1);
r=0.9; gamma=1.4; bcl=1.; bcr=1.;
fprintf(' Example 6.5  \n'); fprintf('  Space step - %f  \n',h); fprintf(' \n');
dd(1:nx-1)=zeros(1,nx-1); ud(1:nx-1)=zeros(1,nx-1); pd(1:nx-1)=zeros(1,nx-1);
x(1:nx)=xmin+((1:nx)-1)*h; xa=xmin+((1:nx-1)-0.5)*h;

%limiter='limiter_mina';
limiter='limiter_supera';
%limiter='limiter_ultraa';
%limiter='limiter_vanAlbada';
%limiter='limiter_vanLeer';
%limiter='limiter_Lin';

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
dd(1:(nx-1)/2)=rol; ud(1:(nx-1)/2)=ul; pd(1:(nx-1)/2)=pl;
dd((nx-1)/2+1:nx-1)=ror;ud((nx-1)/2+1:nx-1)=ur; pd((nx-1)/2+1:nx-1)=pr;

[d,u,p,e,time]=waf_tvd_gasdynamics(dd,ud,pd,h,nx,r,tp,gamma,bcl,bcr,limiter);

% Exact solution
for n=1:nx
    [density, velocity, pressure]=riemann(x(n),time,rol,ul,pl,ror,ur,pr,gamma);
    de(n)=density; ue(n)=velocity; pe(n)=pressure;
    ee(n)=pressure/((gamma-1.)*density);
end
figure(1);
plot(x,de,'k-',xa,d,'k--');
xlabel(' x ');
ylabel(' density ');
title([' Example 6.5 ',te]);
figure(2);
plot(x,ue,'k-',xa,u,'k--');
xlabel(' x ');
ylabel(' velocity ');
title([' Example 6.5 ',te]);
figure(3);
plot(x,pe,'k-',xa,p,'k--');
xlabel(' x ');
ylabel(' pressure ');
title([' Example 6.5 ',te]);
figure(4);
plot(x,ee,'k-',xa,e,'k--');
xlabel(' x ');
ylabel(' internal energy ');
title([' Example 6.5 ',te]);
clear all;