
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    fourth-order finite difference method for Poisson Equation
%-----------------------------------------------------------------------
%    This program solves Poisson equation
%      -( u_{xx} + u_{yy} ) = f(x,y),  a < x <b,  c < y < d.
%    with Dirichlet boundary condition along x=a, x=b, y=c, y=d.
%    
%    Modified by Hanquan Wang. May 25, 2014

clear;  close all

a = 0; b=2.0; c = 0; d=2.0;

m=8*4; n=m;

hx = (b-a)/m; hx1 = hx*hx; x=zeros(m+1,1);
for i=1:m+1,
  x(i) = a + (i-1)*hx;
end

hy = (d-c)/n; hy1 = hy*hy; y=zeros(n+1,1);
for i=1:n+1,
  y(i) = c + (i-1)*hy;
end

%-------B2 -----A2-------------------------------------------------------------
A2 = zeros(m-1,m-1);
B2 = zeros(m-1,m-1) ;

  for i=2:m-2
 A2(i,i) = 10;A2(i,i-1) = 1; A2(i,i+1)= 1;
  end
  
  A2(1,1) =10;
  A2(1,2) = 1;
  
  A2(m-1,m-2)= 1;
  A2(m-1,m-1)= 10;
  
  
for i=2:m-2
  B2(i,i) = -2; B2(i,i-1) = 1; B2(i,i+1)= 1;
end


B2(1,1) = -2; 
B2(1,2) = 1; 

B2(m-1,m-1) = -2; 
B2(m-1,m-2) = 1;

  %-------------------------------------------------------------
  mm=(m-1)*(n-1);
  AA2=sparse(mm,mm);
  BB2=sparse(mm,mm);
  
  I=eye(n-1,n-1);
  
  AA2 =kron(I,A2);
  BB2= kron(I,B2);
  
   
  %-------B1 -----A1-------------------------------------------------------------

B1 = zeros(n-1,n-1) ;
A1 = zeros(n-1,n-1);

  for i=2:n-2
 A1(i,i) = 10;A1(i,i-1) = 1; A1(i,i+1)= 1;
  end
  
  A1(1,1) =10;
  A1(1,2) = 1;
  
  A1(n-1,n-2)= 1;
  A1(n-1,n-1)= 10;
  
  
for i=2:n-2
  B1(i,i) = -2; B1(i,i-1) = 1; B1(i,i+1)= 1;
end


B1(1,1) = -2; 
B1(1,2) = 1; 

B1(n-1,n-1) = -2; 
B1(n-1,n-2) = 1;

  %-------------------------------------------------------------
  mm=(m-1)*(n-1);
 
  AA1=sparse(mm,mm);
  BB1=sparse(mm,mm);
  
  I1=eye(m-1,m-1);
  
  AA1 =kron(A1,I1);
  
  BB1= kron(B1,I1);
  
  M = (m-1)*(n-1); AAA = sparse(M,M); F = zeros(M,1); uu= zeros(M,1);
  %---------------------
  
  AAA = 12.0/(hx^2)*inv(AA2)*BB2 +  12.0/(hy^2)*inv(AA1)*BB1;
  
  %---------------------
  
 
  for j = 2: n,
   for i = 2: m,
    k = i-1 + (j-2)*(m-1);
    F(k) = f(x(i),y(j)); % The right handside of poisson equation
   end
  end
  
  uu = -inv(AAA)*F; % the approximated solutions

  %---------------------
  
  

u= zeros(m+1,n+1);
u2= zeros(m+1,n+1);

% The boundary conditions are dealt with here
for i=1:n+1
  u(1,i) = sin(pi*x(1))*sin(pi*y(i));
  u(m+1,i) =  sin(pi*x(m+1))*sin(pi*y(i));
   u2(1,i) =  sin(pi*x(1))*sin(pi*y(i)); % the exact solution
  u2(m+1,i) =  sin(pi*x(m+1))*sin(pi*y(i)); % the exact solution
end
  
for  i=1:m+1
  u(i,1) =  sin(pi*x(i))*sin(pi*y(1));
  u(i,n+1) =  sin(pi*x(i))*sin(pi*y(n+1));
  u2(i,1) =  sin(pi*x(i))*sin(pi*y(1)); % the exact solution
  u2(i,n+1) =  sin(pi*x(i))*sin(pi*y(n+1)); % the exact solution
end



j = 2;
for k=1:M
  i = k +1 - (j-2)*(m-1) ;
  u(i,j) = uu(k);
  u2(i,j) = ue(x(i),y(j));
  j = fix(k/(m-1)) + 2;
end



% Analyze abd Visualize the result.

e = max( max( abs(u-u2)))        % The maximum error
x1=x(1:m+1); y1=y(1:n+1);

mesh(y1,x1,u); title('近似解'); xlabel('y'); ylabel('x');
figure(2); mesh(y1,x1,u-u2); title('近似解与准确解之差'); xlabel('y');
ylabel('x');
