
%%%%%%%% Clear all unwanted variable and graphs.

clear;  close all

%%%%%%% Input 

a =0; b=1; n=40;
ua = 0; ub = 0;

%%%%%% Call solver: U is 

[x,U] = two_point(a,b,ua,ub,'f',n);

%%%%%%%%%%%%%%%%%%% Plot and error analysis: %%%%%%%%%%%%%%%%%%%
figure(1);
plot(x,U,'o'); % draw the first picture
                        % The approximate solution 

u=zeros(n-1,1);
for i=1:n-1,
  u(i) = sin(pi*x(i));
%  u(i) = x(i)*x(i);
end
hold on;plot(x,u)             % The exact solution

%%%%%%% Plot error

figure(2); plot(x,U-u)

%norm(U-u,inf)		%%% Print out the maximum error.
