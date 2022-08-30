% Example 4.3

l=1.; nx=11; ny=nx;
kappa=1.; ac=1.; tp=0.3;
c=10.; f0=1.;
uu(1:nx,1:ny)=zeros(nx,ny); ud(1:nx,1:ny)=zeros(nx,ny);
hx=l/(nx-1); hy=l/(ny-1);
for n=1:nx
    x(n)=(n-1)*hx;
end
for m=1:ny
    y(m)=(m-1)*hy;
end
tau=0.5*min(hx,hy)^2/kappa;
bx=kappa*tau/(hx*hx); bxa=2.*bx; by=kappa*tau/(hy*hy); bya=2.*by;
fprintf(' Example 4.3  \n'); fprintf('   \n');
fprintf(' Explicit scheme \n'); fprintf('  Space steps - %f   %f  \n',hx,hy);
fprintf('  Time step  - %f  \n',tau);
nts=round(tp/tau); t(1)=0.; ua(1)=0.;
for k=1:nts
    time=k*tau; t(k+1)=time;
    [uu]=heat_2d_es(ud,kappa,ac,nx,ny,hx,hy,tau,time,1,1,'g1_e43','g3_e43',...
                              1,1,'g2_e43','g4_e43','fs_e43');
    ud(1:nx,1:ny)=uu(1:nx,1:ny);
    ua(k+1)=uu((nx-1)/2,(ny-1)/2);
end

% Exact solution
for k=1:nts+1
    time=(k-1)*tau;
    sum=0.;
    for n=0:10
        for m=0:10
            a=pi*n; a2=a*a; b=pi*m; b2=b*b;
            gamma=kappa*(a2+b2)/(l*l); ut=(exp(-gamma*time)-exp(-c*time))/(c-gamma);
            if n == 0
                if m == 0
                    gnm=4.*f0/9.;
                else
                    gnm=-16.*f0*(cos(b)+1.)/(3.*b2);
                end    
            else
                if m == 0
                    gnm=-16.*f0*(cos(a)+1.)/(3.*a2);
                else
                    gnm=64.*f0*(cos(a)+1.)*(cos(b)+1.)/(a2*b2);
                end                    
            end    
            sum=sum+gnm*ut*cos(0.5*a)*cos(0.5*b);
        end
    end
    ue(k)=sum;
end
plot(t,ue,'k-',t,ua,'k--');
xlabel(' t ');
ylabel(' u ');
title(' Example 4.3 ');
clear all;