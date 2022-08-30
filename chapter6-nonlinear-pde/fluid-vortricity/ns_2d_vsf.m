% Two-dimensional Navier-Stokes equation
% Vorticity-Stream function formulation

function [u,v,w,psi,time]=ns_2d_vsf(u,v,w,Re,tp,lx,ly,nx,ny,bxl,bxr,byl,byr,g1,g2,g3,g4)

hx=lx/(nx-1); hy=ly/(ny-1); hxa=1./hx; hya=1./hy; ta=0.5*Re*min(hx,hy)^2;
s(1:nx)=zeros(nx,1); as(1:nx)=zeros(nx,1); bs(1:nx)=zeros(nx,1);
cs(1:nx)=zeros(nx,1);
q(1:ny)=zeros(ny,1); aq(1:ny)=zeros(ny,1); bq(1:ny)=zeros(ny,1);
cq(1:ny)=zeros(ny,1);
x(1:nx)=((1:nx)-1.)*hx; y(1:ny)=((1:ny)-1.)*hy;
time=0.;
while time <= tp
     maxu=max(max(abs(u))); maxv=max(max(abs(v)));
     muv=max(maxu,maxv);
     if muv > 0.1
         tb=1.6*min(hx,hy)/muv; tau=min(ta,tb);
     else
         tb=min(hx,hy); tau=min(ta,tb);
     end
    time=time+tau; fprintf(' time - %f   \n',time);
    
% Runge-Kutta scheme of order 3
    for m=2:ny
        s(:)=w(:,m);
        for n=2:nx-1
            as(n)=u(n,m)*0.5*hxa+hxa^2/Re;
            bs(n)=-2.*hxa^2/Re;
            cs(n)=-u(n,m)*0.5*hxa+hxa^2/Re;
        end
        ds1=tdmv(as,bs,cs,s); ds2=tdmv(as,bs,cs,s+0.5*tau*ds1);
        ds3=tdmv(as,bs,cs,s+0.75*tau*ds2);
        s=s+(2.*ds1+3.*ds2+4.*ds3)*tau/9.;
        w(2:nx-1,m)=s(2:nx-1)';
    end
    for n=2:nx
        q(:)=w(n,:);
        for m=2:ny-1
            aq(m)=v(n,m)*0.5*hxa+hxa^2/Re;
            bq(m)=-2.*hxa^2/Re;
            cq(m)=-v(n,m)*0.5*hxa+hxa^2/Re;
        end
        dq1=tdmv(aq,bq,cq,q); dq2=tdmv(aq,bq,cq,q+0.5*tau*dq1);
        dq3=tdmv(aq,bq,cq,q+0.75*tau*dq2);
        q=q+(2.*dq1+3.*dq2+4.*dq3)*tau/9.;
        w(n,2:ny-1)=q(2:ny-1);
    end
    [psi]=poisson_ns_fft(w,nx,ny,hx,hy);
    for m=2:ny-1
        u(1,m)=0.; u(nx,m)=0.;
        if bxl == 0
            w(1,m)  =-2.*hxa^2*(psi(2,m)+hx*feval(g1,y(m),time));
            v(1,m)=feval(g1,y(m),time);
        else
            w(1,m)=0.; v(1,m)=v(2,m);
        end
        if bxr == 0
            w(nx,m)=-2.*hxa^2*(psi(nx-1,m)-hx*feval(g3,y(m),time));
            v(nx,m)=feval(g3,y(m),time);
        else
            w(1,m)=0.; v(nx,m)=v(nx-1,m);
        end        
    end    
    for n=2:nx-1
        v(n,1)=0.; v(n,ny)=0.;
        if byl == 0
            w(n,1)  =-2.*hya^2*(psi(n,2)-hy*feval(g2,x(n),time));
            u(n,1)=feval(g2,x(n),time);
        else
            w(n,1)=0.; u(n,1)=u(n,2);
        end
        if byr == 0
            w(n,ny)=-2.*hya^2*(psi(n,ny-1)+hy*feval(g4,x(n),time));
            u(n,ny)=feval(g4,x(n),time);
        else
            w(n,ny)=0.; u(n,ny)=u(n,ny-1);
        end
    end    
    for n=2:nx-1    
        for m=2:ny-1        
            u(n,m)= 0.5*(psi(n,m+1)-psi(n,m-1))*hxa;
            v(n,m)=-0.5*(psi(n+1,m)-psi(n-1,m))*hya;
        end
    end
    
   
  figure;  mesh(x,y,psi); drawnow
    figure;  mesh(x,y,u); drawnow
      figure;  mesh(x,y,v); drawnow
      pause;
      
    
  
    
end 
return;