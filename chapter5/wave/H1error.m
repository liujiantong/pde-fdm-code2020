function [L2err,H1err]=H1error(coord,h,uh,u,udx,gamma,T)
nvert=max(size(coord)); 
x=[]; 
k=0;
for i = 1:nvert-1   
   xm = (coord(i+1)+coord(i))*0.5;   
   x = [x, coord(i),xm];   
   k = k + 2;
end
ndof = k+1; 
x (ndof) = coord (nvert); 
uq = eval(u);
uxq = eval(udx);
L2err = 0; 
H1err = 0;
for i=1:nvert-1
   L2err = L2err + (h/6)*((uh(i)-uq(2*i-1))^2+...
           4*(0.5*uh(i)+0.5*uh(i+1)-uq(2*i))^2+(uh(i+1)-uq(2*i+1))^2);
   H1err = H1err + (1/(6*h))*((uh(i+1)-uh(i)-h*uxq(2*i-1))^2+...
           4*(uh(i+1)-uh(i)-h*uxq(2*i))^2+(uh(i+1)-uh(i)-h*uxq(2*i+1))^2);
end
H1err = sqrt(H1err+L2err); 
L2err = sqrt(L2err);
return