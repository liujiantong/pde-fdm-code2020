%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    An Introduction to Scientific Computing          %%%%%%%
%%%%%%%    I. Danaila, P. Joly, S. M. Kaber & M. Postel     %%%%%%%
%%%%%%%                 Springer, 2005                      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Computes the explicit terms for the            %
%  passive scalar equation                        %
%     Hc_s=-(d(u*sca)/dx + d(v*sca)/dy)           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hcs=NSE_calc_hcs(u,v,sca)

global dx dy
global im ip jp jm ic jc

hcs = -0.5/dx*(u(ip,jc).*(sca(ip,jc)+sca) - u.*(sca(im,jc)+sca))...
      -0.5/dy*(v(ic,jp).*(sca(ic,jp)+sca) - v.*(sca(ic,jm)+sca));