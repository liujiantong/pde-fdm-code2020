
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

A2 = zeros(n-1,n-1) ;
M2 = zeros(n-1,n-1);

A = zeros(n-1,n-1) ;
F = zeros(n-1,1);
        
for i=2:n-2
 % A(i,i) = -2/h1; A(i+1,i) = 1/h1; A(i,i+1)= 1/h1;
  A2(i,i) = -24; A2(i,i-1) = 12; A2(i,i+1)= 12;
end

 alpha=-12;
% 
aa = (11*alpha+35)/12;
bb = -(5*alpha+26)/3;
cc = (alpha+19)/2;
dd = (alpha-14)/3;
ee= (11-alpha)/12;

A2(1,1) = aa; % 145/12;
A2(1,2) = bb; %-76/12;
A2(1,3) = cc; %29/2;
A2(1,4) = dd; %-4/3;
A2(1,5) = ee; % 1/12;


A2(n-1,n-1) = aa; %145/12;
A2(n-1,n-2) = bb; %- 76/12;
A2(n-1,n-3) = cc; %29/2;
A2(n-1,n-4) = dd; %-4/3;
A2(n-1,n-5) = ee; % 1/12;

% alpha=11;
% A2(1,1) = 13;
% A2(1,2) = -27;
% A2(1,3) = 15;
% A2(1,4) = -1;
% 
% % 
% % 
% A2(n-1,n-1) = 13;
% A2(n-1,n-2) = - 27;
% A2(n-1,n-3) = 15;
% A2(n-1,n-4) = -1;

% alpha=-5/2;
% ee = -1/4;
% A2(1,1) = alpha+2-ee;
% A2(1,2) = -(2*alpha+5+4*ee);
% A2(1,3) = alpha+4+6*ee;
% A2(1,4) = -(1+4*ee);
% 
% % 
% % 
% A2(n-1,n-1) = alpha+2-ee;
% A2(n-1,n-2) =  -(2*alpha+5+4*ee);
% A2(n-1,n-3) = alpha+4+6*ee;
% A2(n-1,n-4) =-(1+4*ee);




cond(A2)


%  A2(n,n) = -2;
  
%    A2(1,n) = 1;
%    A2(n,1) = 1;
  
  for i=2:n-2
 % A(i,i) = -2/h1; A(i+1,i) = 1/h1; A(i,i+1)= 1/h1;
  M2(i,i) = 10; M2(i,i-1) = 1; M2(i,i+1)= 1;
  end
  
  M2(1,1) =1;
  M2(1,2) =alpha;
  
  M2(n-1,n-2)=alpha;
  M2(n-1,n-1)=1;
  
 cond(M2)
  
A = -1.0/h1*( inv(M2)*A2 );
  
%cond(A)
  
  
  
for i=1:n-1,
  x(i) = a+i*h;
  F(i) = feval(f,x(i)); % 确定函数f(x)在x(i)处的值
end

%U = -h1*( inv(A2)*M2)*F;
  
U = A\F; % 找到方程组的根

return

