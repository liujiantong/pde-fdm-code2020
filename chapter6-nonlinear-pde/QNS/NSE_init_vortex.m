%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    An Introduction to Scientific Computing          %%%%%%%
%%%%%%%    I. Danaila, P. Joly, S. M. Kaber & M. Postel     %%%%%%%
%%%%%%%                 Springer, 2005                      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Initial condition for the vortex dipole       %
% computes the velocity field for a single vortex %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                        % first vortex      
% [rhs,phi]=NSE_init_vortex(Lx,Ly,xc,yc,Lx/4,Ly/2+0.05,0.,0,1);
% u=u+rhs;v=v+phi;
%                        % second vortex
% [rhs,phi]=NSE_init_vortex(Lx,Ly,xc,yc,Lx/4,Ly/2-0.05,-0.,0,-1);
% u=u+rhs;v=v+phi;

function [u,v]=NSE_init_vortex(Lx,Ly,x,y,xv,yv,uin,vin,pm)

global nxm nym dx dy
global im ip jp jm ic jc

% xv, yv  : coordinates of the vortex center
% uin,vin : initial translation velocity

            % radius of the vortex
lv=min([xv,Lx-xv,yv,Ly-yv])*0.400*sqrt(2.);
            % intensity $\psi_0$
psi0=0.1;
            % stream-function and velocity field
for jy=1:nym
for jx=1:nxm
  uloc=(x(jx)-xv)^2+(y(jy)-yv)^2;
  uloc=pm*psi0*exp(-uloc/lv^2);
  u(jx,jy)=uin-2.d0*(y(jy)-yv)*uloc/lv^2;
  v(jx,jy)=vin+2.d0*(x(jx)-xv)*uloc/lv^2;
end
end
