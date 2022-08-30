% Simulating the 1-D Burgers' equation by the Finite Difference
...Method (a time march)
% Numerical scheme used is a first order 
...upwind in time and space for the convection terms 
...and a second order central difference in space
...for the diffusive term

%%
%Specifying Parameters
nx=20;              %Number of steps in space(x)
nt=50;              %Number of time steps 
dt=0.01;            %Width of each time step
dx=2/(nx-1);        %Width of space step
x=0:dx:2;           %Range of x (0,2) and specifying the grid points
u=zeros(1,nx);      %Preallocating u
un=zeros(1,nx);     %Preallocating un
vis=0.1;            %Diffusion coefficient/viscosity
ip=zeros(1,nx);     %Auxillary variable
im=zeros(1,nx);     ...same as above
phi=zeros(1,nx);    ...same as above
dphi=zeros(1,nx);   ...same as above

%%
%Setting up auxillary variables 
for i=1:nx
    ip(i)=i+1;
    im(i)=i-1;
    phi(i)=exp(-0.25*(x(i)^2)/vis)+exp(-0.25*(((2*pi)-x(i))^2)/vis);
    dphi(i)=(-0.5*x(i)/vis)*exp(-0.25*(x(i)^2)/vis)+(0.5*((2*pi)-x(i))/vis)*exp(-0.25*(((2*pi)-x(i))^2)/vis);
end
ip(nx)=1;
im(1)=nx;

t=zeros(nt+1,1);

for i=1:nt+1
    t(i) = 0+(i-1)*dt;
end
%%
%Initial conditions
for i=1:nx
    u(i)=(-2*vis*(dphi(i)/phi(i)))+4;
end

uu= zeros(nx, nt+1);

uu(:,1) = u;
%%
%Explicit scheme with F.D in time and B.D and C.D in space
for it=2:nt+1
    un=u;
    h=plot(x,u);       %plotting the velocity profile
    axis([0 2 4 6])
    title({['1-D Burgers'' equation (\nu = ',num2str(vis),')'];['time(\itt) = ',num2str(dt*it)]})
    xlabel('Spatial co-ordinate (x) \rightarrow')
    ylabel('Transport property profile (u) \rightarrow')
    drawnow;
    refreshdata(h)
    
    for i=1:nx
        if(un(i)>0) 
          u(i)=un(i)-(un(i)*dt*(un(i)-un(im(i)))/dx)+(vis*dt*(un(ip(i))-2*un(i)+un(im(i)))/(dx*dx));
        else
            u(i)=un(i)-(un(i)*dt*(un(ip(i))-un(i))/dx)+(vis*dt*(un(ip(i))-2*un(i)+un(im(i)))/(dx*dx)); 
        end
%         if (x(i)==2*pi)   %boundary condition of periodicity
%             u(x(i))=u(1);
%         end
    end
    
    uu(:,it)=u;
end

figure;mesh(t(3:nt+1),x,uu(:,3:nt+1));