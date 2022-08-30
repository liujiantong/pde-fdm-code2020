clear; close all
 a = 0;  b=1; tfinal = 0.2
 m = 40; 

 h = (b-a)/m; k = 0.8*h; mu = k/h;

 t = 0; n = fix(tfinal/k);
 y1 = zeros(m+1,1); y2=y1; x=y1;   

 for i=1:m+1,
   x(i) = a + (i-1)*h;         % ½Úµã
   y1(i) = u0(x(i));           % ³õÖµ
   y2(i) = 0;
 end

figure(1); plot(x,y1); hold; axis([0 1 -0.2 1.2]);
pause(2)