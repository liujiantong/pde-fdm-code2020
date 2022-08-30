
function [x,U] = two_point(a,b,ua,ub,f,n)

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

A = sparse(n-1,n-1);
F = zeros(n-1,1);
        
for i=1:n-2,
 % A(i,i) = -2/h1; A(i+1,i) = 1/h1; A(i,i+1)= 1/h1;
  A(i,i) = 2/h1; A(i+1,i) = -1/h1; A(i,i+1)= -1/h1;
end
  A(n-1,n-1) = 2/h1;

for i=1:n-1,
  x(i) = a+i*h;
  F(i) = feval(f,x(i)); % 确定函数f(x)在x(i)处的值
end
  F(1) = F(1)+ua/(2.0*h);
  F(n-1) = F(n-1) +ub/(2.0*h);
	
U = A\F; % 找到方程组的根

return

