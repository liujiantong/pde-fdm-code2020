function advectiondiffusion
% File name: advectiondiffusion.m
% Author: Clive Mingham 
% Date: 14 Oct 2010
%
% This m-file is associated with the free text book: 'Introductory Finite Difference 
% Methods for PDEs' which may be downloaded from
% http://bookboon.com/uk/student/mathematics/introductory-finite-difference
% -methods-for-pdes
%
% Description: Solves the 1D Advection-Diffusion equation 
% dU/dt + v dU/dx = Kx d2U/dx2 with conditions given by Figure 5.2.
%
% Variables:
% u=u(x,t) is the concentration function
% runtime is run time in seconds
% t is elapsed time in seconds
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
% Sub function: gaussian
%
% Boundary Conditions:
% Dirichlet - set to ghost values to zero at left and right edges of domain
%
% Output:
% The solution is calculated over [p, q] using N points and plotted
% after ntimesteps time steps, i.e. the final 
% solution is at time, ntimesteps*dt seconds.
%
clear all; clc; clf
% First set the initial parameters.
p=0;                           % start of computational domain
q=100;                         % end of computational domain
v=0.5;                         % water speed in x direction
t=0;                           % start clock
runtime=57;                    % required run time
Kx=0.1;                        % diffusion coeficient in x direction
N=101;                         % number of grid points
dx=(q-p)/(N-1);                % spatial step size
x = p: dx : q;                 % vector of grid points
u0=zeros(1,N);                 % vector of initial u values filled with zeros
% Set initial profile
u0=gaussian(x);                 % gaussian initial profile
initialu=u0;                    % keep initial profile for comparison
%
F=0.9;                         % safety factor
dtAD=dx*dx/(v*dx+2*Kx);        % unsplit time step from von Neumann stability analysis
dt=F*dtAD;                     % time step reduced by safety factor
ntimesteps=100;                % number of time steps
Cx=v*dt/dx;                   % Courant number
Rx=Kx*(dt/(dx*dx));            % constant for diffusion term
% 
u=zeros(1,N+2);                % define correct sized numerical solution array
%                              
% Begin the time marching algorithm
disp('start time marching ...')
 for count=1:ntimesteps
     t=t+dt;        % current time for outputted solution
     if t > runtime
         t=t-dt;    % go back
         dt=t-runtime; % calc last dt
         t=runtime
     else
         t
     end                                 
    u0=[0 u0  0];         % insert ghost values
    for i=2:N+1
        u(i)=u0(i)-Cx*(u0(i)-u0(i-1))+ Rx*(u0(i+1)-2*u0(i)+u0(i-1));   
    end
    u0=u(2:N+1);            % copy solution to initial conditions for next iteration
    if t >= runtime
         disp('runtime achieved')
         break   % stop time marching loop     
     end      
   end  % of time marching loop 
   disp('time marching finished, see results in graphics window')
   % Output
    plot(x,initialu,'k:',x,u(2:N+1),'k+')  % plot of numerical soln and initial profile
    xlabel('x')
    ylabel('concentration u')
    title('advection-diffusion: initial profile (dotted) and unsplit numerical solution (+) at later time, ') 
% end of advectiondiffusion.sce
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = gaussian(x)
  % Initial condition Gaussian function over [0,100], height h, centred on x = xcen
  N=length(x);
  y=x*0;
  h=1.0;
  xcen=20;
  for i=1:N
      if (x(i)>=10 && x(i)<=30) 
        y(i)=h*exp(-0.05*(x(i)-xcen)^2); 
      end    
  end  
% end of function gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


