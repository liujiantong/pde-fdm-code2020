
clear all; close all
 a = 0;  b=1; tfinal = 0.5
 m = 200; 

 h = (b-a)/m; k = 0.2*h; mu = 2.0*k/h;

 t = 0; n = fix(tfinal/k);
 y1 = zeros(m+1,1); y2=y1; x=y1;   

 
 for i=1:m+1,
   x(i) = a + (i-1)*h;         % 节点
   y1(i) = u0(x(i));           % 初值
   y2(i) = 0;
 end

%figure(1); plot(x,y1); hold;
%pause(2)

%u=zeros(n,m+1);
 

 t = 0;        % Starting time loop.
 for j=1:n,
   y1(1)=uexact(t,x(1)); y2(1)=uexact(t+k,x(1));      % boundary conditions
   for i=2:m+1
     y2(i) = y1(i) - mu*(y1(i)-y1(i-1) );
      
   end

   %plot(x,y2); pause(0.5)
   
   t = t + k; y1 = y2;
 end

 u_e = zeros(m+1,1);
 for i=1:m+1
   u_e(i) = uexact(t,x(i));
 end

 max(abs(u_e-y2))

 figure(2); plot(x,y2,':',x,u_e,'-b');
 xlabel('x'); title('t=0.5时的近似解与准确解');
figure(3); plot(x,u_e-y2);
    xlabel('x'); title('t=0.5时近似解与准确解之差');

