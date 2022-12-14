
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    
%    This program solves Poisson equation
%       u_{xx} + u_{yy} = f(x,y),  a < x <b,  c < y < d.
%    with zero Dirichlet boundary condition along x=a, x=b, y=c, y=d.
%       ???߽?????
%    Modified by Hanquan Wang. Oct. 25 2009

% clear;  close all

function u=poisson_matlab(m,n,f)

a = 0; b=15; c = -15; d=15;

%m=150; n=150;

hx = (b-a)/m; hx1 = hx*hx; x=zeros(m+1,1);
for i=1:m+1,
  x(i) = a + (i-1)*hx;
end

hy = (d-c)/n; hy1 = hy*hy; y=zeros(n+1,1);
for i=1:n+1,
  y(i) = c + (i-1)*hy;
end


M = (m-1)*(n-1); A = sparse(M,M); bf = zeros(M,1);
% We use k = i-1 + (j-2)*(m-1), i=2,3,...m; j=2,3,...,n
% to index the program.

%f=zeros(m+1,n+1);

for j = 2: n,
  for i = 2: m,
    k = i-1 + (j-2)*(m-1);
    %bf(k) = f(x(i),y(j)); % The right handside of poisson equation
    bf(k) = f(i,j);
    A(k,k) = -2/hx1 -2/hy1;

    %-- x direction --------------

    if i == 2
        A(k+1,k) = 1/hx1;
     %   bf(k) = bf(k) - ue(a,y(j))/hx1;
    else
       if i== m
         A(k-1,k) = 1/hx1;
      %   bf(k) = bf(k) - ue(b,y(j))/hx1;
       else
          A(k-1,k) = 1/hx1;
          A(k+1,k) = 1/hx1;
       end
     end

%-- y direction --------------

    if j == 2
        A(k,k+m-1) = 1/hy1;
       % bf(k) = bf(k) - ue(x(i),c)/hy1;
    else
       if j== n
         A(k,k-m+1) = 1/hy1;
        % bf(k) = bf(k) - ue(x(i),d)/hy1;
       else
          A(k,k-m+1) = 1/hy1;
          A(k,k+m-1) = 1/hy1;
       end
     end

  end
end

  U = A \bf;

%--- Transform back to (i,j) form to plot the solution ---

u= zeros(m+1,n+1);
%u2= zeros(m+1,n+1);

% % The boundary conditions are dealt with here
% for i=1:n+1
%   u(1,i) = exp(x(1))*sin(pi*y(i));
%   u(m+1,i) = exp(x(m+1))*sin(pi*y(i));
%    u2(1,i) = exp(x(1))*sin(pi*y(i)); % the exact solution
%   u2(m+1,i) = exp(x(m+1))*sin(pi*y(i)); % the exact solution
% end
%   
% for  i=1:m+1
%   u(i,1) = exp(x(i))*sin(pi*y(1));
%   u(i,n+1) = exp(x(i))*sin(pi*y(n+1));
%   u2(i,1) = exp(x(i))*sin(pi*y(1)); % the exact solution
%   u2(i,n+1) = exp(x(i))*sin(pi*y(n+1)); % the exact solution
% end

j = 2;
for k=1:M
  i = k +1 - (j-2)*(m-1) ;
  u(i,j) = U(k);
 % u2(i,j) = ue(x(i),y(j));
  j = fix(k/(m-1)) + 2;
end

% Analyze abd Visualize the result.

%e = max( max( abs(u-u2)))        % The maximum error
%x1=x(1:m+1); y1=y(1:n+1);

%mesh(y1,x1,u); title('The solution plot'); xlabel('y'); ylabel('x');
%figure(2); mesh(y1,x1,u-u2); title('The error plot'); xlabel('y');
%ylabel('x');

end