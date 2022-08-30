
%%%%%%%% Clear all unwanted variable and graphs.

clear;  close all

%%%%%%% Input 

a =1; b=2; N=50;

h = (b-a)/N;

x = zeros(N+1,1);
y = zeros(N+1,1);

for n=1:N+1,
  x(n) = a+(n-1)*h;
end

y(1) =0; % ��ʼ����

% y_{n+1} = y_n + h*f(x_n,y_n), n=0,1,2....N-1 % ��ǰEuler����ʽ


for n=1:N
   
    y(n+1) = y(n) + h*f(x(n),y(n)); % f��x(n),y(n))ȷ������f(x,y)��x(n),y(n)����ֵ
    
end


figure(1);
plot(x,y,'o'); % draw the first picture
                        % The approximate solution 

y_exact=zeros(N+1,1);
for n=1:N+1,
%  u(i) = cos(pi*x(i));
  y_exact(n) = x(n)*x(n)*(exp(x(n))-exp(1.0));
end
hold on;plot(x,y_exact)             % The exact solution

%%%%%%% Plot error

figure(2); plot(x,y_exact-y)

norm(y_exact-y,inf)		%%% Print out the maximum error.
