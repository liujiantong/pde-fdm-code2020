%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    An Introduction to Scientific Computing          %%%%%%%
%%%%%%%    I. Danaila, P. Joly, S. M. Kaber & M. Postel     %%%%%%%
%%%%%%%                 Springer, 2005                      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Computes and prints the divergence of the     %
%   velocity field                                %
%     div(q)=du/dx+dv/dy                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NSE_print_div(u,v,niter,temps)

global dx dy
global im ip jp jm ic jc
global ncarl ;

      divns = (u(ip,jc)-u)/dx+ (v(ic,jp)-v)/dy;

divmax=max(max(abs(divns)));    % max of divergence
divtot=sqrt(sum(sum(divns.*divns))*dx*dy); % integral

fprintf('It=%d   time=%5.3f  div_max=%10.5e div_tot=%10.5e  \n',niter,temps,divmax,divtot);