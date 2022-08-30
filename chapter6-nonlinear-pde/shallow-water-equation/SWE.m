function SWE
% File name: SWE.m
% Author Clive Mingham 
% Date: 19 Oct 2010
%
% This m-file is associated with the free text book: 'Introductory Finite Difference 
% Methods for PDEs' which may be downloaded from
% http://bookboon.com/uk/student/mathematics/introductory-finite-difference
% -methods-for-pdes
%
% Description: Solves the 1D Shallow Water (SWE) equation, dU/dt + dF/dx = 0, where,
% U = [h    ,  F = [      h*v
%      h*v]          h*v^2 + 0.5*g*h^2] 
% using the Lax-Friedrichs scheme. 
% The solution, U, is calculated over [p, q] using N points and h and/or v are plotted.
% Set up parameters are for the column collapse in Chapter 7.
% Define variables and their [units]:
% runtime = the required run time [s]
% t = current time [s]
% dt = current time step [s]
% x = distance along the flow [m]
% safe = safety factor for time step
% h = 1 by N+2 matrix for height [m] including ghost values
% v = 1 by N+2 matrix for velocity [m/s] incuding ghost values
% U0 = 2 by N+2 matrix for solution at time level n including ghost values
% U = 2 by N+2 matrix for solution at time level n+1 including ghost values 
% F = 2 by N+2 matrix for fluxes at time level n including ghost values
%
% Define constants:
% g = the acceleration due to gravity [m/s^2]
%
% Boundary Conditions: Need left and right ghost values for Lax-Friedrichs scheme so put variables in 
% arrays whose columns go from i=1 to i=N+2.  
% Left and right ghost indices are i=1 and i=N+2 respectively. Computational values go from i=2 to i=N+1.
% LBC = left boundary condition flag (1 for solid, 0 for transmissive)
% RBC = right boundary condition flag (1 for solid, 0 for transmissive)
%
% Note: This program does not reproduce the original Figure 7.1 which was
% generated by an incorrect program.  
%
% This program has been verified on a small test problem but it would be a worthwhile
% exercise to verify it on your own small problem.
%
% Subfunctions: initialh, initialv,fillU, fillF, calcdt, gethandv 

clear all; clc; clf
% Set parameters:
g=9.81;  
safe=0.95; % safety factor for time step
p=0;       % start of computational domain
q=100;     % end of computational domain
LBC=0;     % flags left boundary condition
RBC=0;     % flags right boundary condition
t=0.0;     % start time
runtime=1.0;  % run time [s]
N=201;         % number of grid points
dx=(q-p)/(N-1);  % spatial step size
x=p:dx:q;        % vector of grid points    
ntimesteps=0;    % time step counter        
U0=zeros(2,N+2);  % define U0 array - time level n
U=U0;             % define U array - time level n+1
% Set initial h and v (in the computational domain)
h=initialh(x); 
hinit=h;          % store for plotting
v=initialv(x);
vinit=v;          % store for plotting
%
disp('program running...') 
% Begin the time marching algorithm
while (t<runtime)  % time loop
     ntimesteps=ntimesteps+1;
     dt=calcdt(h,v,dx);        % calculate timestep
     dt=safe*dt;               % use safety factor
     t=t+dt;                   % current time
%   Implement boundary conditions for h and v
%   Left boundary
     if (LBC==1)    % solid left boundary
        h(1)=h(2); 
        v(1)=-v(3);    
     elseif (LBC==0)  % transmissive left boundary
            h(1)=h(2);
            v(1)=v(2);
     else
           disp('no left boundary conditions so stop')
           exit
     end  % of left boundary conditions           
%   Right boundary
     if (RBC==1)     % solid right boundary
        h(N+2)=h(N+1);
        v(N+2)=-v(N);
     elseif (RBC==0)  % transmissive right boundary
            h(N+2)=h(N+1);
            v(N+2)=v(N+1);
     else
           disp('no right boundary conditions so stop')
           exit
         end  % of right boundary conditions  
     % Fill U0 and F values which are 2x(N+2) matrices
     U0=fillU(h,v);
     F=fillF(h,v);    
%   Implement FD scheme - computational grid points start and i=2 and end at i=N+1
     for i=2:N+1
         U(1,i)=0.5*(U0(1,i+1)+U0(1,i-1))-(dt/(2*dx))*(F(1,i+1)-F(1,i-1));
         U(2,i)=0.5*(U0(2,i+1)+U0(2,i-1))-(dt/(2*dx))*(F(2,i+1)-F(2,i-1));
      end  
      [h,v]=gethandv(U);   % extract h and v at computational grid points for next iteration      
end  % of time loop
disp('program finished, see graphics window')
% Output
subplot(2,1,1)
plot(x,hinit(2:N+1),'k--',x,h(2:N+1))  % plot of height and initial height (--)
title('Lax-Friedrichs numerical solution to the Shallow Water Equations')
xlabel('x [m]')
ylabel('h water depth [m]')
%
subplot(2,1,2)
plot(x,vinit(2:N+1),'k--',x,v(2:N+1))  % plot of velocity and initial velocity (--)
title('Lax-Friedrichs numerical solution to the Shallow Water Equations')
xlabel('x [m]')
ylabel('v water velocity [m/s]')
%
disp('number of time step =')
ntimesteps
disp('run time=')
t
% end of program SWE.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = initialh(x)
% Initial height function
% models a dam break into a lake in 1D 
  N=length(x);    % number of computational grid points
  h=zeros(1,N+2); % append left (i=1) and right (i=N+2) ghost values to be determined from boundary conditions
                  % and set to zero for now.
  depthlake=1.0;    % depth of lake outside column [m] 
  depthcol=10.0;    % depth of water column      
  xstartcol=40.0;   % x coordinate of start of water column [m]
  xendcol=60.0;     % x coordinate of end of water column [m]                          
  for i=1:N
      if (x(i)>=xstartcol & x(i)<=xendcol) 
         h(i+1)=depthcol;  % note computational indices start at i=2
       else
         h(i+1)=depthlake;
      end 
  end 
% end of function initialh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = initialv(x)
% Initial velocity function   
  N=length(x);
  v=zeros(1,N+2);  % append left and right ghost values to be determined from boundary conditions
% end of function initialv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = fillU(h,v)
% function to fill U - include ghost values
  U(1,:)=h;
  U(2,:)=h.*v;
% end of function fillU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = fillF(h,v)
% function to fill F - include ghost values
  g=9.81;
  F(1,:)=h.*v;
  F(2,:)=h.*v.*v+0.5*g*h.*h;
% end of function fillF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dt = calcdt(h,v,dx)
% function to find stable time step dt
  g=9.81;
  celerity=sqrt(g*h);
  vel1=max(abs(v+celerity));
  vel2=max(abs(v-celerity));
  maxwavespeed=max(vel1,vel2);
  dt=dx/maxwavespeed;
% end of function calcdt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h,v] = gethandv(U)
% function to get h and v from U at the computational grid points i=2 to i=N+1
% ghost values are are implemented later
% computational values go from i=2
  drythreshold=0.01;  % below this [m] the bed is considered dry
  M=length(U(1,:));  % number of colums in U including ghost entries
  N=M-2;  % number of computational grid points
  h=zeros(1,M);  % create h with space for left and right ghost values
  v=zeros(1,M);  % create v with space for left and right ghost values
  h(2:N+1)=U(1,2:N+1);
  % loop needed for v in case h = 0
  for i=2:N+1
    if (h(i)<drythreshold)  % no water, avoids division by zero
       v(i)=0;
    else
      v(i)=U(2,i)/U(1,i); 
    end
  end
% end of function gethandv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
