% Simulating the 2-D Linear Convection equation by the Finite Difference
...Method 
% Numerical scheme used is a first order upwind in both space and time

%%
%Specifying Parameters
nx=90;                           %Number of steps in space(x)
ny=90;                           %Number of steps in space(y)       
nt=50;                           %Number of time steps 
dt=0.01;                         %Width of each time step
c=1.0;                           %Velocity of wave propagation
dx=2/(nx-1);                     %Width of space step(x)
dy=2/(ny-1);                     %Width of space step(y)
x=0:dx:2;                        %Range of x(0,2) and specifying the grid points
y=0:dy:2;                        %Range of y(0,2) and specifying the grid points
u=zeros(nx,ny);                  %Preallocating u
un=zeros(nx,ny);                 %Preallocating un

%%
%Initial Conditions
for i=1:nx
    for j=1:ny
        if ((0.5<=y(j))&&(y(j)<=1)&&(0.5<=x(i))&&(x(i)<=1))
            u(i,j)=2;
        else
            u(i,j)=0;
        end
    end
end

%%
%Boundary conditions
u(1,:)=0;      
u(nx,:)=0;
u(:,1)=0;
u(:,ny)=0;

%%
%Explicit method with F.D in time and B.D in space
i=2:nx-1;
j=2:ny-1;
for it=1:nt
    un=u;
    h=surf(x,y,u','EdgeColor','none');       %plotting the velocity profile
    shading interp
    axis([0 2 0 2 0 2.5])
    title(['2-D Linear Convection with {\itc} = ',num2str(c)])
    xlabel('Spatial co-ordinate (x) \rightarrow')
    ylabel('{\leftarrow} Spatial co-ordinate (y)')
    zlabel('Transport property profile (u) \rightarrow')
    drawnow;
    refreshdata(h)
    u(i,j)=un(i,j)-(c*dt*(un(i,j)-un(i-1,j))/dx)-(c*dt*(un(i,j)-un(i,j-1))/dy);
    %Boundary conditions
    u(1,:)=0;
    u(nx,:)=0;
    u(:,1)=0;
    u(:,ny)=0;
end