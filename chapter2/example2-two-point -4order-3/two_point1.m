
function [x,U] = two_point1(a,b,ua,ub,f,n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    This matlab function two_point solvers the following two-point    %
%    boundary value problem: -u''(x) = f(x) using the center difference %
%    scheme.                                                           %  
%                                                                      %
%    Input:                                                            %
%     a, b: Two end points.                                            %
%     ua, ub: Dirichlet boundary conditions at a and b                 %
%     f: external function f(x).                                       %
%     n: number of grid points.                                        %
%                                                                      %
%    Output:                                                           %
%     x: x(1),x(2),...x(n-1) are grid points                           %
%     U: U(1),U(2),...U(n-1) are approximate solution at the grid      %
%     points   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = (b-a)/n; h1=h*h;

A2 = zeros(n-1,n-1) ;
M2 = zeros(n-1,n-1);
temp = zeros(n-1,1);

A = zeros(n-1,n-1) ;
F = zeros(n-1,1);
        
for i=2:n-2
 % A(i,i) = -2/h1; A(i+1,i) = 1/h1; A(i,i+1)= 1/h1;
  A2(i,i) = -2; A2(i,i-1) = 1; A2(i,i+1)= 1;
end


A2(1,1) = -2; % 145/12;
A2(1,2) = 1; %-76/12;

A2(n-1,n-1) = -2; % 145/12;
A2(n-1,n-2) = 1; %-76/12;

cond(A2)
  
  for i=2:n-2
 % A(i,i) = -2/h1; A(i+1,i) = 1/h1; A(i,i+1)= 1/h1;
  M2(i,i) = 10; M2(i,i-1) = 1; M2(i,i+1)= 1;
  end
  
  M2(1,1) =10;
  M2(1,2) = 1;
  
  M2(n-1,n-2)= 1;
  M2(n-1,n-1)= 10;

cond(M2)
  
   
for i=1:n-1,
  x(i) = a+i*h;
  F(i) = feval(f,x(i)); % 确定函数f(x)在x(i)处的值
end

temp = -M2*F;

temp(1) = temp(1) -(12/h1+1.0)*feval(f,a);
temp(n-1) = temp(n-1) -(12/h1+1.0)*feval(f,b);


temp = h1/12*temp;

U = inv(A2)*temp; % 找到方程组的根

return

