
%%%%%%%% Clear all unwanted variable and graphs.

clear;  close all

%%%%%%% Input 

a =-1; b= 1; n=40;
ua = 0; ub = 0;

h = (b-a)/n; h1=h*h;

x = zeros(n+1,1); u=x; 
v=zeros(n+1,1);
       
for i=1:n+1,
  x(i) = a+(i-1)*h;
end

eps=1.0e-6;
error =1.0;
while (error>eps)
   for i=2:n
    u(i) = 0.5*( -h1*exp(v(i))+(v(i+1)+v(i-1))); 
   end 
   u(1)=0; 
   u(n+1)=0;
   error =max(abs(u-v));
   v=u;
   error
end




%%%%%%%%%%%%%%%%%%% Plot and error analysis: %%%%%%%%%%%%%%%%%%%
figure(1);
plot(x,u,'-o'); hold on; % draw the first picture
                        % The approximate solution 

