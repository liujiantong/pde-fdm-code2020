
%%%%%%%% Clear all unwanted variable and graphs.

clear;  close all

%%%%%%% Input 

a =-1; b= 1;c=-1;d=1;
n=40;
ua = 0; ub = 0;
uc=0; ud =0;

hx = (b-a)/n; hx1=hx*hx;

hy = (d-c)/n; hy1=hy*hy;




x = zeros(n+1,1); 
y=x;

u=zeros(n+1,n+1);
v=zeros(n+1,n+1);
       
for i=1:n+1,
  x(i) = a+(i-1)*hx;
end

for j=1:n+1,
  y(j) = c+(j-1)*hy;
end


eps=1.0e-6;
error =1.0;
while (error>eps)
   for i=2:n
       for j=2:n
    u(i,j) = 0.25*(-hx1*exp(v(i,j))+(v(i+1,j)+v(i-1,j))+...
         (v(i,j+1)+v(i,j-1)));
     
%  u(i,j) = (-hx1*exp(v(i,j))+(v(i+1,j)+v(i-1,j))+...
%          (v(i,j+1)+v(i,j-1))*hx1/hy1)/(2.0+2.0*hx1/hy1);
     
       end 
   end
   
   error =max(max(abs(u-v)));
   v=u;
   error
end




%%%%%%%%%%%%%%%%%%% Plot and error analysis: %%%%%%%%%%%%%%%%%%%
figure(1);
[X,Y] = meshgrid(x', y');
mesh(X,Y, u);


