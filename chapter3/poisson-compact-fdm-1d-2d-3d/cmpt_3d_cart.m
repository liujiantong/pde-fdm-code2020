% to test the compact finite difference scheme for 3D Laplace equation
% (u_xx + u_yy + u_zz) = f, (x,y,z) \in [a,b]^3;
%
% test two examples 1) u = sin(pi/(b-a) (x-a)) sin( pi/(b-a)(y-a))sin(
% pi/(b-a)(z-a))
%			  2) u = exp( - (x^2+y^2+z^2) )
%
N = 2^2;  a = -8; b= 8; x = linspace(a,b,N+1);x = x(1:end-1); y = x;z=x;
h= (b-a)/N; mu_x = pi*(1:(N-1))/(b-a);  mu_x= mu_x'; mu_y = mu_x; mu_z=mu_x;
tmp = 2*(cos(mu_x*h)-1)/h^2;  mu_mul = tmp./(1+h^2/12*tmp);

f = zeros(N,N,N); u = f;
for i =1:N
	for j = 1:N
      for k = 1:N
	u_ext(i,j,k) = sin(pi/(b-a)*(x(i)-a))*sin(pi/(b-a)*(y(j)-a))...
                  *sin(pi/(b-a)*(z(k)-a));
	f(i,j,k) = -u_ext(i,j,k)*3.0*(pi/(b-a))^2;

%	u_ext(i,j) =exp(-(x(i)^2+y(j)^2+z(k)^2));
%	f(i,j) = 4*u_ext(i,j)*(x(i)^2+y(j)^2-1);
    
      end
	end
end

ftmp = f(2:N,2:N,2:N);
ftmp1 =ftmp;

%f_dst= dst(dst(dst(ftmp)')')'; 

for i =1:N-1
	for j = 1:N-1
        ftmp1(i,j,:)= dst(ftmp1(i,j,:));
    end 
end

for i =1:N-1
	for j = 1:N-1
        ftmp(i,:,j)= dst(ftmp1(i,:,j));
    end 
end

for i =1:N-1
	for j = 1:N-1
        ftmp1(:,i,j)= dst(ftmp(:,i,j));
    end 
end


u_dst = ftmp;
for k =1:N-1
	for l = 1:N-1
        for p = 1:N-1
	u_dst(k,l,p) = ftmp1(k,l,p)/(mu_mul(k)+ mu_mul(l)+mu_mul(p));
        end
    end
end

for i =1:N-1
	for j = 1:N-1
        ftmp1(i,j,:)= idst(u_dst(i,j,:));
    end 
end

for i =1:N-1
	for j = 1:N-1
        ftmp(i,:,j)= idst(ftmp1(i,:,j));
    end 
end

for i =1:N-1
	for j = 1:N-1
        ftmp1(:,i,j)= idst(ftmp(:,i,j));
    end 
end

u(2:N,2:N,2:N) =ftmp1;
%u(2:N,2:N,2:N) =idst(idst(idst(u_dst)')')';

max(max(max(abs(u_ext(2:N,2:N,2:N)-u(2:N,2:N,2:N)))))

