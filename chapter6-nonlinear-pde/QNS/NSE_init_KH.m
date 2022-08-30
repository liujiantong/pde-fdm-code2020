%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    An Introduction to Scientific Computing          %%%%%%%
%%%%%%%    I. Danaila, P. Joly, S. M. Kaber & M. Postel     %%%%%%%
%%%%%%%                 Springer, 2005                      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Initial condition for a 2D periodic jet       %
%    study of the Kelvin-helmholtz  instability   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u=NSE_init_KH(Lx,Ly,x,y,U0,Rj,Pj,Ax,lamx)

global dx dy
global im ip jp jm ic jc
           %  2D grid
[xx,yy]=meshgrid(x,y);xx=xx';yy=yy';
           % local $y$ coordinate
rr=yy-Ly/2;
           % velocity $u$
u=U0*0.5d0*(1.d0+tanh(0.5*Pj*(1-abs(rr)/Rj)));
           % perturbation in $x$ direction
u=u+Ax*u.*sin(2*pi/lamx*xx);

% sca=NSE_init_KH(Lx,Ly,xm,ym,1,Ly/4,20,0.00,0.5*Lx);