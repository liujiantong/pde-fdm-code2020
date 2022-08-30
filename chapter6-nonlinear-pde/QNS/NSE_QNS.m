%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    An Introduction to Scientific Computing          %%%%%%%
%%%%%%%    I. Danaila, P. Joly, S. M. Kaber & M. Postel     %%%%%%%
%%%%%%%                 Springer, 2005                      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Resolution of 2D Navier-Stokes equations       %
%   incompressible fluid                          %
%=================================================%
%  rectangular 2D domain (L_x,L_y)                % 
%  with periodic boundary conditions in  x and y  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     close all; clear all;
     format long e;
%===============================%
%  Global variables             %
%===============================%
     global dx dy Lx Ly;
     global nxm nym ;
     global ip im jp jm ic jc;

%===============================%
%  Input parameters             %
%  Initial condition            %
%===============================%
icas=input(['请从1、2、3、4之中选择运行例子\n'... 
                            '1=Kelvin-Helmholtz (1)\n'...
                            '2=Kelvin-Helmholtz (2)\n'...
                            '3=dipole (1)\n'...
                            '4=dipole (2)\n']);
                        
% icas=input(['Choose the run case\n'... 
%                             '1=Kelvin-Helmholtz (1)\n'...
%                             '2=Kelvin-Helmholtz (2)\n'...
%                             '3=dipole (1)\n'...
%                             '4=dipole (2)\n']);
                                                

                        
switch icas
case {1,2}
  %fprintf('Evolution of the Kelvin-Helholtz instability\n');  
  fprintf('Kelvin-Helholtz 非稳定性的运行过程\n'); 
      Lx=2; Ly=1;
      nx=65;ny=65;
      rey=1000;
      pec=1000;

      nitermax=210;
      nprint=10;niso=10;  
   case{3,4}
   % fprintf('Evolution of a vortex dipole \n');  
    fprintf('涡旋的运动过程 \n'); 
      Lx=1; Ly=1;
      nx=65;ny=65;
      rey=1000;
      pec=1000;

      nitermax=500;
      nprint=100;niso=10;
   end  

%===============================%
%  2D grid                      %
%===============================%
      nxm=nx-1  ;           nym=ny-1;
      dx=Lx/nxm ;           dy=Ly/nym;

      ic=1:nxm;             jc=1:nym; 

      xc=(ic-1)*dx ;        yc=(jc-1)*dy ;
      xm=(ic-0.5)*dx;       ym=(jc-0.5)*dy;

      ip=ic+1; ip(nxm)=1;   jp=jc+1; jp(nym)=1;
      im=ic-1; im(1)=nxm;   jm=jc-1; jm(1)=nym;
                             
[xx,yy]=meshgrid(xm,ym);xx=xx';yy=yy'; % centers of the cells for visualization

%===============================%
%  Initialization               %
%===============================%

      u   =zeros(nxm,nym);      % velocity u
      v   =zeros(nxm,nym);      % velocity v

      gpu  =zeros(nxm,nym);     % pressure gradient along x
      gpv  =zeros(nxm,nym);     % pressure gradient along y
      hcu  =zeros(nxm,nym);     % explicit terms for u
      hcv  =zeros(nxm,nym);     % explicit terms for v
      hcs  =zeros(nxm,nym);     % explicit terms for the passive scalar

      pres =zeros(nxm,nym);     % pressure
      sca  =zeros(nxm,nym);     % passive scalar
     
      rhs  =zeros(nxm,nym);     % work array for the RHS
      phi  =zeros(nxm,nym);     % variable for the pressure correction (Poisson eq)

%===============================%
%  Initial condition            %
%===============================%
switch icas
case 1
                       % initial velocity field
      u  =NSE_init_KH(Lx,Ly,xc,ym,1,Ly/4,20,0.25,0.5*Lx);
      sca=NSE_init_KH(Lx,Ly,xm,ym,1,Ly/4,20,0.00,0.5*Lx);
                       % cfl for the time step 
      cfl=0.2;
case 2
                       % initial velocity field
      u  =NSE_init_KH(Lx,Ly,xc,ym,1,Ly/4,20,0.25,0.25*Lx);
      sca=NSE_init_KH(Lx,Ly,xm,ym,1,Ly/4,20,0.00,0.25*Lx);
                       % cfl for the time step 
      cfl=0.1;
case 3
    
                           % one vortex      
[rhs,phi]=NSE_init_vortex(Lx,Ly,xc,yc,Lx/4,Ly/2,0.,0,1);
u=u+rhs;v=v+phi;
%                        % first vortex of one pair     
% [rhs,phi]=NSE_init_vortex(Lx,Ly,xc,yc,Lx/4,Ly/2+0.05,0.,0,1);
% u=u+rhs;v=v+phi;
%                        % second vortex of one pair
% [rhs,phi]=NSE_init_vortex(Lx,Ly,xc,yc,Lx/4,Ly/2-0.05,-0.,0,-1);
% u=u+rhs;v=v+phi;
                       % scalar field (a stripe in the middle of the domain)
      sca(nxm/2-10:nxm/2+10,:)=ones(21,nym);
                       % cfl for the time step 
      cfl=0.4;
case 4
                        % first pair of vortices -> dipole 1
[rhs,phi]=NSE_init_vortex(Lx,Ly,xc,yc,Lx/4,Ly/2+0.05,0.,0,1);
u=u+rhs;v=v+phi;
[rhs,phi]=NSE_init_vortex(Lx,Ly,xc,yc,Lx/4,Ly/2-0.05,-0.,0,-1);
u=u+rhs;v=v+phi;
                        % second pair of vortices -> dipole 2
[rhs,phi]=NSE_init_vortex(Lx,Ly,xc,yc,3*Lx/4,Ly/2+0.05,-0.,0,-1);
u=u+rhs;v=v+phi;
[rhs,phi]=NSE_init_vortex(Lx,Ly,xc,yc,3*Lx/4,Ly/2-0.05,-0.,0,1);
u=u+rhs;v=v+phi;
                       % scalar field (a stripe in the middle of the domain)
      sca(nxm/2-10:nxm/2+10,:)=ones(21,nym);
                       % cfl for the time step
      cfl=0.4;        
end  

% visualization of the initial condition
%=========================================%
ivisu=1;
figure(1);
%set(1,'MenuBar','none');
     NSE_visu_vort(xc,yc,u,v,11,ivisu,0);
% figpos=get(1,'Position');
% figure(2);set(2,'MenuBar','none');
% set(2,'Position',[figpos(1) figpos(2)-figpos(4)-20 figpos(3) figpos(4)]);
% NSE_visu_sca (xm,ym,sca,11,ivisu,0);

pause 
%===============================%
% Optimization of the ADI method%
%===============================%
      dt=NSE_calc_dt(u,v,cfl);       % time step
    %  fprintf('Time step dt=%d \n',dt);  
       fprintf('时间步长 dt=%d \n',dt);
      
      bx=0.5*dt/(dx*dx)/rey;
      by=0.5*dt/(dy*dy)/rey;

[amix,apix,alphx,xs2x]=NSE_ADI_init(-bx*ones(1,nxm),(1+2*bx)*ones(1,nxm),-bx*ones(1,nxm));
[amiy,apiy,alphy,xs2y]=NSE_ADI_init(-by*ones(1,nym),(1+2*by)*ones(1,nym),-by*ones(1,nym));

%=========================================%
%   Optimization of the Poisson solver    %
%=========================================%

      kl=(cos(2*pi/nxm*(ic-1))-1)*2/(dx*dx);

      am=ones(nxm,nym)/(dy*dy);
      ac=(-2/(dy*dy)+kl)'*ones(1,nym);
      ap=ones(nxm,nym)/(dy*dy);

          % wave-number = 0
      am(1,1)=0;
      ac(1,1)=1;
      ap(1,1)=0;
     
      [am,ap,ac,xs2]=NSE_Phi_init(am,ac,ap);
      
%===============================%
%   Time loop                   %
%===============================%

niter=0; temps=0;
NSE_print_div(u,v,niter,temps);

tc=cputime; % to estimate the computational CPU time: initialization
   while(niter <= nitermax)
      
         niter=niter+1;temps=temps+dt;

      % momentum equation for  u
      %===========================
                                     % RHS term
	 rhs =-0.5*hcu;
	 hcu = NSE_calc_hcu(u,v);
	 rhs = dt*(rhs+1.5*hcu-gpu+NSE_calc_lap(u)/rey);
                                     % first step of ADI
	 du1 = NSE_ADI_step(amix,apix,alphx,xs2x,rhs');
                                     % second step of ADI
	 du  = NSE_ADI_step(amiy,apiy,alphy,xs2y,du1');
                                     % non-solenoidal u field
          u = u+du;

      % momentum equation for  v
      %===========================
                                     % RHS term
	 rhs =-0.5*dt*hcv;
	 hcv = NSE_calc_hcv(u,v);
	 rhs = dt*(rhs+1.5*hcv-gpv+NSE_calc_lap(v)/rey);
                                     % first step of ADI
	 du1 = NSE_ADI_step(amix,apix,alphx,xs2x,rhs');
                                     % second step of ADI
	  du = NSE_ADI_step(amiy,apiy,alphy,xs2y,du1');
                                     % non-solenoidal v field
          v=v+du;

      % compute the divergence of the non-solenoidal field
      %===========================
  
          rhs=((u(ip,jc)-u)/dx+ (v(ic,jp)-v)/dy)/dt;

      % resolution of the Poisson equation
      %===========================

          uf=fft(rhs); uf(1,1)=0;
          uff=NSE_Phi_step(am,ap,ac,xs2,uf);
          phi=real(ifft(uff));

      % velocity field correction     
      %===========================    

          u=u-dt*(phi-phi(im,jc))/dx;
          v=v-dt*(phi-phi(ic,jm))/dy;

      % compute the pressure field    
      %===========================    

          pres=pres+phi-dt/(2*rey)*NSE_calc_lap(phi);

      % update the pressure gradient    
      %===========================    

          gpu=(pres-pres(im,jc))/dx;
          gpv=(pres-pres(ic,jm))/dy;

      % solve the passive scalar equation
      %=================================
                                     % RHS term
	 rhs =-0.5*hcs;
	 hcs = NSE_calc_hcs(u,v,sca);
	 rhs = dt*(rhs+1.5*hcs+NSE_calc_lap(sca)/pec);
                                     % first step of ADI
	 du1 = NSE_ADI_step(amix,apix,alphx,xs2x,rhs');
                                     % second step of ADI
	  du = NSE_ADI_step(amiy,apiy,alphy,xs2y,du1');
                                     % scalar field
          sca=sca+du;

      % visualization of the flow (vorticity+scalar)
      %==============================

         if(mod(niter,nprint) == 0); 
             NSE_print_div(u,v,niter,temps);ivisu=ivisu+1;
             figure(1);NSE_visu_vort(xc,yc,u,v,niso,ivisu,temps);
             %fname=['print -depsc figV' num2str(ivisu)];eval(fname);
            % figure(2);NSE_visu_sca (xm,ym,sca,niso,ivisu,temps);
             pause 
             %fname=['print -depsc figS' num2str(ivisu)];eval(fname);
          end;

    end;
fprintf('CPU time =%d\n',cputime-tc);
