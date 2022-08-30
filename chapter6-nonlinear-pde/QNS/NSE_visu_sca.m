%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    An Introduction to Scientific Computing          %%%%%%%
%%%%%%%    I. Danaila, P. Joly, S. M. Kaber & M. Postel     %%%%%%%
%%%%%%%                 Springer, 2005                      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Visualization of iso-contours                 %
%   of the passive scalar (or tracer)             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NSE_visu_sca(x,y,sca,ncont,ivisu,temps)

global dx dy
global im ip jp jm ic jc


[xx,yy]=meshgrid(x,y);xx=xx';yy=yy';
pcolor(xx,yy,sca);axis equal;shading('interp');drawnow
set(gca,'FontSize',16);xlabel('x');ylabel('y');
title(['Scalar 函数 在 t=' num2str(temps),'时的值'],'FontSize',16);