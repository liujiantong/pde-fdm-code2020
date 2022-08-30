%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    An Introduction to Scientific Computing          %%%%%%%
%%%%%%%    I. Danaila, P. Joly, S. M. Kaber & M. Postel     %%%%%%%
%%%%%%%                 Springer, 2005                      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Visualization of iso-contours                 %
%   of the vorticity                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NSE_visu_vort(x,y,u,v,ncont,ivisu,temps)
%NSE_visu_vort(xc,yc,u,v,niso,ivisu,temps);

global dx dy
global im ip jp jm ic jc
            % compute the vorticity
omeg=(v(ic,jc)-v(im,jc))/dx-(u(ic,jc)-u(ic,jm))/dy ;
            % 2D grid
[xx,yy]=meshgrid(x,y);xx=xx';yy=yy';
            % iso-contours
           % pcolor(xx,yy,omeg);axis equal;shading('interp');drawnow
          % figure; surfl(xx,yy,omeg);shading flat;colormap gray;
          %  figure; mesh(xx,yy,u);
          %   figure; mesh(xx,yy,v);
             figure;mesh(xx,yy,omeg);
            
            
            
%set(gca,'FontSize',16);
xlabel('x');ylabel('y');
%title(['Vorticity 函数在 t=' num2str(temps),'时的值'],'FontSize',16);