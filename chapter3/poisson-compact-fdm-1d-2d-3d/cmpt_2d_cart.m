% to test the compact finite difference scheme for 2D Laplace equation
% (u_xx + u_yy ) = f, (x,y) \in [a,b]^2;
%
% test two examples 1) u = sin(pi/(b-a) (x-a)) sin( pi/(b-a)(y-a))
%			  2) u = exp( - (x^2+y^2) )
%
N = 2^5;  a = -8; b= 8; x = linspace(a,b,N+1);x = x(1:end-1); y = x;
h= (b-a)/N; mu_x = pi*(1:(N-1))/(b-a);  mu_x= mu_x'; mu_y = mu_x; 
tmp = 2*(cos(mu_x*h)-1)/h^2;  mu_mul = tmp./(1+h^2/12*tmp);

f = zeros(N,N); u = f;
for i =1:N
	for j = 1:N

	u_ext(i,j) = sin(pi/(b-a)*(x(i)-a))*sin(pi/(b-a)*(y(j)-a));
	f(i,j) = -u_ext(i,j)*2*(pi/(b-a))^2;

%	u_ext(i,j) =exp(-(x(i)^2+y(j)^2));
%	f(i,j) = 4*u_ext(i,j)*(x(i)^2+y(j)^2-1);

	end
end

ftmp = f(2:N,2:N);
f_dst= dst(dst(ftmp)')'; u_dst = f_dst;
for k =1:N-1
	for l = 1:N-1
	u_dst(k,l) = f_dst(k,l)/(mu_mul(k)+ mu_mul(l));
	end
end
u(2:N,2:N) = idst(idst(u_dst)')';
max(max(abs(u_ext(2:N,2:N)-u(2:N,2:N))))

