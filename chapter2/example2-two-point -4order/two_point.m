
function [x,U] = two_point(a,b,ua,ub,f,n)


h = (b-a)/n; h1=h*h;

A2 = sparse(n-1,n-1) ;
M3 = sparse(n-1,n-1);

A = sparse(n-1,n-1) ;
F = zeros(n-1,1);
        
for i=1:n-2,
 % A(i,i) = -2/h1; A(i+1,i) = 1/h1; A(i,i+1)= 1/h1;
  A2(i,i) = -2; A2(i+1,i) = 1; A2(i,i+1)= 1;
end
  A2(n-1,n-1) = -2;
  
%   A2(1,n-1) = 1;
%   A2(n-1,1) = 1;
  
  for i=1:n-2,
 % A(i,i) = -2/h1; A(i+1,i) = 1/h1; A(i,i+1)= 1/h1;
  M3(i,i) = 10; M3(i+1,i) = 1; M3(i,i+1)= 1;
  end
   M3(n-1,n-1) = 10;

A = -12/h1* ( (M3)^(-1)*A2 );
  
  
for i=1:n-1
  x(i) = a+i*h;
  F(i) = feval(f,x(i)); % 确定函数f(x)在x(i)处的值
end
  
U = A\F; % 找到方程组的根

return

