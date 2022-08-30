%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    An Introduction to Scientific Computing          %%%%%%%
%%%%%%%    I. Danaila, P. Joly, S. M. Kaber & M. Postel     %%%%%%%
%%%%%%%                 Springer, 2005                      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Computes the explicit terms for the            %
%  v component of velocity                        %
%      Hc_v=-(d(uv)/dx + d(v^2)/dy)               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hcv=NSE_calc_hcv(u,v)

global dx dy
global im ip jp jm ic jc

hcv = -0.25/dx*((u(ip,jc)+u(ip,jm)).*(v+v(ip,jc)) - (u(ic,jm)+u).*(v(im,jc)+v))...
      -0.25/dy*((v(ic,jp)+v).*(v(ic,jp)+v)-(v(ic,jm)+v).*(v(ic,jm)+v));
