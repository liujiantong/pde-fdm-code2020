
%%%%%%%% Clear all unwanted variable and graphs.

clear;  close all

%%%%%%% Input 

a =1; b=2; N=40*1;

h = (b-a)/N;

x = zeros(N+1,1);
y = zeros(N+1,1);

for n=1:N+1,
  x(n) = a+(n-1)*h;
end

y(1) =0; % 初始条件

% y_{n+1} = y_n + h*f(x_n,y_n), n=0,1,2....N-1 % 向前Euler法公式


for n=1:N
   
    %y(n+1) = y(n) + h*f(x(n),y(n)); %
    %f（x(n),y(n))确定函数f(x,y)在x(n),y(n)处的值(the forward Euler method)
     K1=f(x(n),y(n));
     K2=f(x(n+1),y(n)+h*K1);          % this method called 2-order Runge-Kutta
     y(n+1)=y(n)+0.5*h*(K1+K2);
end

% for n=1:N
%    y(n+1)=y(n)+h*f(x(n),y(n)); 
% end
% error=1;
% y1=zeros(N+1,1);y1(1)=0;
% while(error>10^(-3))
%     
% for n=1:N                             % this method called Crank-Nicolson
%     y1(n+1)=y(n)+0.5*h*(f(x(n),y(n))+f(x(n+1),y(n+1)));
%     
% end
% error=abs(max(y-y1));
% y=y1;
% 
% end


figure(1);
plot(x,y,'o'); % draw the first picture
                        % The approximate solution 

y_exact=zeros(N+1,1);
for n=1:N+1,
%  u(i) = cos(pi*x(i));
  y_exact(n) = x(n)*x(n)*(exp(x(n))-exp(1.0));
end
hold on;plot(x,y_exact)             % The exact solution
xlabel('x');
%ylabel('近似解与准确解');
legend('近似解','准确解');

%%%%%%% Plot error

figure(2); plot(x,y_exact-y)
xlabel('x');
ylabel('|y_exact-y|');

norm(y_exact-y,inf)		%%% Print out the maximum error.
