%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    An Introduction to Scientific Computing          %%%%%%%
%%%%%%%    I. Danaila, P. Joly, S. M. Kaber & M. Postel     %%%%%%%
%%%%%%%                 Springer, 2005                      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Initial condition for the vortex dipole       %
% Lamb-Chaplygin vortex dipole (see Batchelor)    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u,v]=NSE_init_dipole(xv,yv,lv,x,y)

global nxm nym dx dy
global im ip jp jm ic jc

%----------- psi_0---------------
       psi0=0.5d0;
%----------- k=3.83/a --------------
       uu0=3.83d0/lv;
%----------  U=-k*psi_0/2*J_0(3.83)--
       u0=0.5d0*uu0*psi0*0.402758809533;
%---------- velocity field
       for j=1:nym
        for i=1:nxm
          uloc=(x(i)-xv)*(x(i)-xv)+(y(j)-yv)*(y(j)-yv);
          uloc=sqrt(uloc);      % local radius
          if(uloc ~= 0)
           cth=(x(i)-xv)/uloc;  % cos(theta)
           sth=(y(j)-yv)/uloc;  % sin(theta)
          else
           cth=0.d0;
           sth=0.d0;
          end
                                % stream-function psi
          if(uloc < lv) 
             rho(i,j)=psi0*sth*besselj(1,uu0*uloc);
          else
             rho(i,j)=0.d0;
          end
        end
     end
     
%-----vitesses---------

      u=(rho(ic,jp)-rho(ic,jc))/dy;    % u= d(psi)/dy
      v=-(rho(ip,jc)-rho(ic,jc))/dx;   % v=-d(psi)/dx
