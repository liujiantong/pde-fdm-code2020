
clear; close all
 a = 0;  b=1; tfinal = 0.5
% m = 20; 
 m = 200;
% m = 80;

 h = (b-a)/m; k = 0.2*h; mu = 2.0* k/h; mu1=mu*mu/2.0;

 t = 0; n = fix(tfinal/k);
 y1 = zeros(m+1,1); y2=y1; x=y1;

 for i=1:m+1,
   x(i) = a + (i-1)*h;          % The grid.
   y1(i) = u0(x(i));            % Initial value
   y2(i) = 0;
 end

%figure(1); plot(x,y1); hold;
%pause(2)

 t = 0;
 for j=1:n,                    % Starting time loop.
   y1(1)=uexact(t,x(1)); y2(1)=uexact(t+k,x(1));
   for i=2:m
%    y2(i) = 0.5*(y1(i+1)+y1(i-1)) - mu*(y1(i+1)-y1(i-1))/2;
     y2(i)= y1(i)-mu*(y1(i+1)-y1(i-1))/2 + mu1*(y1(i+1)-2*y1(i)+y1(i-1));
   end
   i = m+1;
   y2(i)=y1(i)-mu*(3*y1(i)-4*y1(i-1)+y1(i-2))/2+mu1*(y1(i)-2*y1(i-1)+y1(i-2));

 %  plot(x,y2); pause(0.5)

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

