function advectiondiffusionsplit
% File name: advectiondiffusionsplit.m  
% Author: Clive Mingham 
% Date: 18 Oct 2010.
%
% This m-file is associated with the free text book: 'Introductory Finite Difference 
% Methods for PDEs' which may be downloaded from
% http://bookboon.com/uk/student/mathematics/introductory-finite-difference
% -methods-for-pdes
%
% Description: Solves the 1D Advection-Diffusion (AD) equation 
% dU/dt + v dU/dx = Kx d2U/dx2.
% by solving:  
% LA:  dU/dt + v dU/dx = 0 
% then 
% LD:  dU/dt = Kx d2U/dx2
% solution is un+2 = LD(dt)LA(dt)un
%
% Variables:
% U=U(x,t) is the concentration function
% runtime is the required run time in seconds
% t is time in seconds
% x is distance in metres
% v is velocity m/s
% Kx is the diffusion coefficient in the x direction
%
% Note:  
% 1) if Kx = 0 this is pure advection and the problem
%    has the usual analytical solution.
% 2) if v = 0 this is pure diffusion and the problem
%    has an analytical solution for certain inital conditions.
% 
% Discretisation:
% first order forward difference for du/dt
% first order backward difference for du/dx
% second order symmetric difference for d2u/dx2
%
% Subfunctions: gaussian
%
% Boundary Conditions:
% Dirichlet - set to ghost values to zero at left and right edges of domain
%
% Output:
% The solution is calculated over [p, q] using N points and plotted
% every time step, dt, for ntimesteps time steps, i.e. the final 
% solution is at time, ntimesteps*dt seconds.
%
clear all 
clc
clf
% First set the initial parameters.
p=0;                           % start of computational domain
q=100;                         % end of computational domain
t=0.0;                         % start time
runtime=57;                    % run time
v=0.5;                         % water speed in x direction
Kx=0.1;                        % diffusion coeficient in x direction
N=101;                         % number of grid points
dx=(q-p)/(N-1);                % spatial step size
x = p: dx : q;                 % vector of grid points
u0=zeros(1,N);                 % vector of initial u values filled with zeros
for i=1:N
    u0(i)=gaussian(x(i));      % gaussian initial profile
    % u0(i)=sin(pi*x(i));        % special IC for pure diffusion
end
initial=u0;                    % keep initial profile for comparison
%
F=0.9;                         % safety factor
dtA=dx/abs(v);                 % pure advection time step
dtD=dx*dx/(2*Kx);              % pure diffusion time step
%dtAD=dx*dx/(v*dx+2*Kx);       % unsplit time step
dt=F*min(dtA,dtD);             % split time step
ntimesteps=100;                % number of time steps
Cx =v*dt/dx;                   % Courant number
Rx=Kx*(dt/(dx*dx));            % constant for diffusion term
% 
u=zeros(1,N+2);                % define correct sized numerical solution array
ubar=u;        % define correct sized numerical solution array for intermediate step
%                              
% Begin the time marching algorithm
%
M = moviein(ntimesteps);       % creates matrix for movie
 for n=1:ntimesteps
     n
     t=t+2*dt;        % current time for outputted solution
     if t > runtime
         t=t-2*dt;    % go back
         dt=t-runtime; % calc last dt
         t=runtime
     else
         t
     end             
%      for i=1:N
%           exact(i)=gaussian(x(i)-v*t);   % exact solution at time t
%           exact(i)=exp(-Kx*pi*pi*t)*sin(pi*x(i));
%      end                      
    u0=[0 u0  0];         % insert ghost values
%   LA step
    for i=2:N+1
        ubar(i)=u0(i)-Cx*(u0(i)-u0(i-1));   
    end
    ubar=[0 ubar  0];         % insert ghost values
%   LD step
    for i=2:N+1
        u(i)=ubar(i)+ Rx*(ubar(i+1)-2*ubar(i)+ubar(i-1));   
    end
% advection-diffusion
    plot(x,initial,'k:',x,u(2:N+1),'k+')  % plot of numerical soln and initial profile
    xlabel('x')
    ylabel('concentration u')
    title('advection-diffusion: initial profile (dotted) and split numerical solution (+) at later time, ')  
%
    pause(0.01)             % pause in seconds between plots
    u0=u(2:N+1);            % copy solution to initial conditions for next iteration 
    if t >= runtime
         disp('runtime achieved')
         break   % stop time loop     
     end         
 end  % of time loop 
end  % of advectiondiffusionsplit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = gaussian(x)
% Gaussian function for initial conditions over [0,100], height h, centred
% on x = xcen 
  h=1.0;
  xcen=20;
  y=0.0;
  if and(x>=10,x<=30) 
     y = h*exp(-0.05*(x-xcen)^2);     
  end  
end  % of function gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
