clear;close all;
N=10;
M=10;
K=100;
a=0;
b=1;
c=0;
d=1;
t=1/2;
f=1/3;
dx=1/N;
dy=1/M;
r=dx^2;
h=dx;
x=0:1/N:1;
y=0:1/M:1;
u=ones(N+1,M+1,K+1);
w=ones(N+1,M+1,K+1);
u(1,:,2)=sin(y); 
u(N+1,:,2)=exp(1)*sin(y);  
u(:,1,2)=0;  
u(:,M+1,2)=exp(x)*sin(1);
%====================???????????????¡§??¡§?===================%
for k=2:K
for j=2:M
    for i=2:N
       u(i,j,k)=(1/4)*(r*w(i,j,k-1)+u(i+1,j,k)+u(i-1,j,k)+u(i,j-1,k)+u(i,j+1,k));
          for m=1:M+1
               for n=1:N+1
                u(n,m,k+1)= t*u(n,m,k+1)+(1-t)*u(n,m,k);
                w(:,:,k)=-u(:,:,k+1);
               end
            end  
    end
end
end
w(1,:,1)=sin(y); 
w(N+1,:,1)=exp(1)*sin(y);  
w(:,1,1)=0;  
w(:,M+1,1)=exp(x)*sin(1);
for k=1:K
 for j=2:M
 for i=2:N
    w(i,j,k+1)=(-1/4)*((h/2)*((u(i,j+1,k+1)-u(i,j-1,k+1))*(w(i+1,j,k)-w(i-1,j,k))-(u(i+1,j,k+1)-u(i-1,j,k+1))*(w(i,j+1,k)-w(i,j-1,k)))-w(i+1,j,k+1)-w(i-1,j,k+1)-w(i,j+1,k+1)-w(i,j-1,k+1));
       
            w(i,j,k+1)= t*w(i,j,k+1)+(1-t)*w(i,j,k);
          
    end
end
end
figure(1)
mesh(u(:,:,K+1))
title('u');
figure(2)
mesh(w(:,:,K))
title('w');
