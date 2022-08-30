
%%%%%%%% Clear all unwanted variable and graphs.

clear;  close all

%%%%%%% Input 

a =1; b=2; N=500;

h = (b-a)/N;

x = zeros(N+1,1);
y = zeros(N+1,1);

for n=1:N+1,
  x(n) = a+(n-1)*h;
end

y(1) =0; % 初始条件

% y_{n+1} = y_n + h*f(x_{n+1},y_{n+1}), n=0,1,2....N-1 % 向后Euler法公式


for n=1:N
     temp = 1.0-h*2/x(n+1);
    y(n+1) =1/temp*( y(n) + h*x(n+1)*x(n+1)*exp(x(n+1))); % % 向后Euler法公式
    
end


figure(1);
plot(x,y,'o'); % draw the first picture
                        % The approximate solution 

y_exact=zeros(N+1,1);
for n=1:N+1,

  y_exact(n) = x(n)*x(n)*(exp(x(n))-exp(1.0));
  
end
hold on;plot(x,y_exact)             % The exact solution

%%%%%%% Plot error

figure(2); plot(x,y_exact-y)

norm(y_exact-y,inf)		%%% Print out the maximum error.
