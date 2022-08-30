%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    An Introduction to Scientific Computing          %%%%%%%
%%%%%%%    I. Danaila, P. Joly, S. M. Kaber & M. Postel     %%%%%%%
%%%%%%%                 Springer, 2005                      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Computes the explicit terms for the            %
%  u component of velocity                        %
%     Hc_u=-(du^2/dx + d(uv)/dy)                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hcu=NSE_calc_hcu(u,v)

global dx dy
global im ip jp jm ic jc

hcu = -0.25/dx*((u(ip,jc)+u).*(u(ip,jc)+u) - (u(im,jc)+u).*(u(im,jc)+u))...
      -0.25/dy*((u(ic,jp)+u).*(v(ic,jp)+v(im,jp))-(u(ic,jm)+u).*(v(im,jc)+v));
 
