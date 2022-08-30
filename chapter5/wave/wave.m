%  Newmark Scheme. 
%Solution and error of the wave equation
%
%    u_tt - c^2*u_xx = 0
%
% x \in I=[a,b]  and  t \in [0,T] with the 
% initial data  u(x,0)  = u_0(x) 
%           and u'(x,0) = v_0(x)
% and the boundary conditions
%   u(a,t) = 0 and u(b,t) = 0 
% with the  Newmark Scheme.
%
% This method is presented for example in 
% Quarteroni/Sacco/Saleri: Numerische Mathematik 2

% Author     : Stefan Hüeber
% Date       : Februar 3, 2003
% Institution: University of Stuttgart,
%              Institut for Applied Analysis and Numerical Mathematics,
%              Numerical Mathematics for super computers  
% Version    : 1.0


clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial data
u0 = inline('2*sin(pi*x)');
% time derivative
v0 = inline('4*pi*gamma*sin(2*pi*x)','x','gamma');
% exact solution
uu = 'sin(pi*(x-gamma*T))+sin(pi*(x+gamma*T))+cos(2*pi*(x-gamma*T))-cos(2*pi*(x+gamma*T))';
ut = ['gamma*pi*(-cos(pi*(x-gamma*T))+cos(pi*(x+gamma*T)))+2*gamma*',...
			' pi*(sin(2*pi*(x-gamma*T))+sin(2*pi*(x+gamma*T)))'];
uu1 = inline('sin(pi*(x-gamma*T))+sin(pi*(x+gamma*T))+cos(2*pi*(x-gamma*T))-cos(2*pi*(x+gamma*T))','x','gamma','T');
ut1 = inline('gamma*pi*(-cos(pi*(x-gamma*T))+cos(pi*(x+gamma*T)))+2*gamma*pi*(sin(2*pi*(x-gamma*T))+sin(2*pi*(x+gamma*T)))','x','gamma','T');

% velocity
gamma = 1;

% Newmark-Scheme
beta = 0.25;
theta = 0.5;

% Intervall 
T = 4;
I = [0 1];

% Number of time steps
TT = 200;
% Number of space steps
XX = 50;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% menu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag = menu('Select',...
						'Newmark',...
						'Errorplot');

disp('Please wait ...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch flag
 case 1
	errsteps = 1;
 case 2
	T = 2;
	% Number of time steps
	TT = 250;
	% Number of space steps
	XX = 2;
	% Newmark-Scheme
	betaopt = 0.25;
	thetaopt = 0.5;
	betanonopt = 0.25;
	thetanonopt = 0.3;
	% Number of refinement steps for error plot
	errstepsopt = 9;
	errstepsnonoptder = 7;
	errstepsnonoptfun = errstepsopt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Euler forward for time / central for space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag == 2
	beta = betaopt;
	theta = thetaopt;
else
	errstepsopt = 1;
end
	

TTopt(1) = TT;
XXopt(1) = XX;
	
for k = 1:errstepsopt
	NT = TTopt(k);
	NX = XXopt(k);
	
	dt = T/NT;
	dx = (I(2)-I(1))/NX;
	lambda = dt/dx;
	
	% Initial conditions
	SOL = zeros(2*(NX+1),2);
	if flag == 1
		SOLPLOT = zeros(NX+1,NT+1);
	end
	for j = 2:NX
		% value of the function
		SOL(j,1) = u0(I(2)+(1-j)*dx);
	end
	if flag == 1
		SOLPLOT(:,1) = SOL(1:NX+1,1);
	end;
	for j = 1:NX+1
		% value of the time derivative
		SOL(NX+1+j,1) = v0(I(2)+(1-j)*dx,gamma);
	end
	
	% Finite Difference-Scheme: Newmark-Scheme
	M = spdiags([ones(NX,1),-2*ones(NX,1),ones(NX,1)],...
							[-1,0,1],NX-1,NX-1);
	MAT = speye(NX-1) - (gamma*lambda)^2*beta*M;
	for n = 1:NT
		% value of the function
		rhs = (speye(NX-1) + (gamma*lambda)^2*(0.5-beta)*M) ...
					* SOL(2:NX,1) + dt*SOL(NX+3:2*NX+1,1);
		SOL(2:NX,2) = MAT\rhs;
		% value of the time derivative
		SOL(NX+3:2*NX+1,2) = ((gamma*lambda)^2/dt) * ...
				(theta*M*SOL(2:NX,2) + (1-theta)*M*SOL(2:NX,1)) + ...
				SOL(NX+3:2*NX+1,1);
		SOL(:,1) = SOL(:,2);
		if flag == 1 
			SOLPLOT(:,n+1) = SOL(1:NX+1,1);
		end
	end
	
	% Computation of the L2-error for t=T
	xx = linspace(I(2),I(1),NX+1);
	[L2erropt(k),H1err(k)] = H1error(xx,dx,SOL(1:NX+1,1),uu,'0*x', ...
																	 gamma,T);
	[L2erroptder(k),H1err(k)] = H1error(xx(2:NX),dx,SOL(NX+3:2*NX+1,1),ut,'0*x', ...
																		 gamma,T);
	
	
	if flag==1
		% Visualization
		time = ones(NX+1,1)*linspace(0,T,NT+1);
		space = linspace(I(2),I(1),NX+1)'*ones(1,NT+1);
		surf(time,space,SOLPLOT(1:NX+1,:),'EdgeColor','none');
		xlabel('Zeit t','FontSize',16);
		ylabel('Ort x','FontSize',16);
		pause(0.5);
	end
	
	if k~=errstepsopt
		XXopt(k+1) = 2*XXopt(k);
		TTopt(k+1) = 2*TTopt(k);
	end
end



if flag == 2	
	beta = betanonopt;
	theta = thetanonopt;
	
	TTnonopt(1) = TT;
	XXnonopt(1) = XX;
	
	% Order of the function for nonoptimal theta
	for k = 1:errstepsnonoptfun
		NT = TTnonopt(k);
		NX = XXnonopt(k);
		
		dt = T/NT;
		dx = (I(2)-I(1))/NX;
		lambda = dt/dx;
		
		% Initial conditions
		SOL = zeros(2*(NX+1),2);
		for j = 2:NX
			% value of the function
			SOL(j,1) = u0(I(2)+(1-j)*dx);
		end
		
		for j = 1:NX+1
			% value of the time derivative
			SOL(NX+1+j,1) = v0(I(2)+(1-j)*dx,gamma);
		end
		
		% Finite Difference-Scheme: Newmark-Scheme
		M = spdiags([ones(NX,1),-2*ones(NX,1),ones(NX,1)],...
								[-1,0,1],NX-1,NX-1);
		MAT = speye(NX-1) - (gamma*lambda)^2*beta*M(1:NX-1,1:NX-1);
		for n = 1:NT
			% value of the function
			rhs = (speye(NX-1) + (gamma*lambda)^2*(0.5-beta)*M(1:NX-1,1:NX-1)) ...
						* SOL(2:NX,1) + dt*SOL(NX+3:2*NX+1,1);
			SOL(2:NX,2) = MAT\rhs;
			% value of the time derivative
			SOL(NX+3:2*NX+1,2) = ((gamma*lambda)^2/dt) * ...
					(theta*M*SOL(2:NX,2) + (1-theta)*M*SOL(2:NX,1)) + ...
					SOL(NX+3:2*NX+1,1);
			SOL(:,1) = SOL(:,2);
		end		
		
		% Computation of the L2-error for t=T
		xx = linspace(I(2),I(1),NX+1);
		[L2errnonopt(k),H1err(k)] = H1error(xx,dx,SOL(1:NX+1,1),uu,'0*x',gamma,T);
		if k~=errstepsnonoptfun
			XXnonopt(k+1) = 2*XXnonopt(k);
			TTnonopt(k+1) = 2*TTnonopt(k);   
		end
	end
	
	% Order of the time derivative for nonoptimal theta
	XXnonoptder(1) = XX;
	TTnonoptder(1) = TT;

	for k = 1:errstepsnonoptder
		NT = TTnonoptder(k);
		NX = XXnonoptder(k);
		
		dt = T/NT;
		dx = (I(2)-I(1))/NX;
		lambda = dt/dx;
		
		% Initial conditions
		SOL = zeros(2*(NX+1),2);
		for j = 2:NX
			% value of the function
			SOL(j,1) = u0(I(2)+(1-j)*dx);
		end
		
		for j = 1:NX+1
			% value of the time derivative
			SOL(NX+1+j,1) = v0(I(2)+(1-j)*dx,gamma);
		end
		
		% Finite Difference-Scheme: Newmark-Scheme
		M = spdiags([ones(NX,1),-2*ones(NX,1),ones(NX,1)],...
								[-1,0,1],NX-1,NX-1);
		MAT = speye(NX-1) - (gamma*lambda)^2*beta*M(1:NX-1,1:NX-1);
		for n = 1:NT
			% value of the function
			rhs = (speye(NX-1) + (gamma*lambda)^2*(0.5-beta)*M(1:NX-1,1:NX-1)) ...
						* SOL(2:NX,1) + dt*SOL(NX+3:2*NX+1,1);
			SOL(2:NX,2) = MAT\rhs;
			% value of the time derivative
			SOL(NX+3:2*NX+1,2) = ((gamma*lambda)^2/dt) * ...
					(theta*M*SOL(2:NX,2) + (1-theta)*M*SOL(2:NX,1)) + ...
					SOL(NX+3:2*NX+1,1);
			SOL(:,1) = SOL(:,2);
		end
		
		% Computation of the L2-error for t=T
		xx = linspace(I(2),I(1),NX+1);
		[L2errnonoptder(k),H1err(k)] = H1error(xx(2:NX),dx,SOL(NX+3:2*NX+1,1),ut,'0*x',...
																					gamma,T);
		if k~=errstepsnonoptder
			XXnonoptder(k+1) = 2*XXnonoptder(k);
			TTnonoptder(k+1) = 4*TTnonoptder(k); 
		end
	end


	figure(1);
	% Errorplot	of function
	loglog(XXopt.*TTopt, 15*((1./TTopt).^2 + (1./XXopt).^2),'c--', ...
				 'linewidth',3);
	hold on; box on;
	loglog(XXnonopt.*TTnonopt,L2errnonopt,'bo-','Linewidth',3, ...
				 'MarkerFaceColor','b');
	loglog(XXopt.*TTopt,L2erropt,'mo--','Linewidth',3, ...
				 'MarkerFaceColor','m');
	title(['L2-error of u for t=',num2str(T)],'FontSize',24,'Color','r');
	xlabel('N_t*N_x','FontSize',16,'Color','k');
	ylabel('Fehler','FontSize',16,'Color','k');
	set(gca,'FontSize',16);
	set(gca,'xlim',[min(XXnonopt.*TTnonopt), max(XXnonopt.* ...
																							 TTnonopt)],'ylim',[0 60]);
	legend(['O(\Delta t^2 + \Delta x^2)'],...
				 ['\theta = ',num2str(thetanonopt),', \Delta t_{n+1} = 0.5\Delta t_n, \Delta x_{n+1} = 0.5\Delta x_n'],...
				 ['\theta = ',num2str(thetaopt),', \Delta t_{n+1} = 0.5\Delta',' t_n, \Delta x_{n+1} = 0.5\Delta x_n'],1);
	
	figure(2);
	% Errorplot	of function
	loglog(XXnonoptder.*TTnonoptder, 15*((1./TTnonoptder) + (1./XXnonoptder).^2),'g--', ...
				 'linewidth',3);
	hold on; box on;
	loglog(XXopt.*TTopt, 15*((1./TTopt).^2 + (1./XXopt).^2),'c--', ...
				 'linewidth',3);
	loglog(XXnonoptder.*TTnonoptder,L2errnonoptder,'mo-','Linewidth',3, ...
				 'MarkerFaceColor','m');
	loglog(XXopt.*TTopt,L2erroptder,'bo-','Linewidth',3, ...
				 'MarkerFaceColor','b');
	title(['L2-error of u_t for t=',num2str(T)],'FontSize',24,'Color','r');
	xlabel('N_t*N_x','FontSize',16,'Color','k');
	ylabel('Fehler','FontSize',16,'Color','k');
	set(gca,'FontSize',16);
	set(gca,'xlim',[min(XXnonopt.*TTnonopt), max(XXnonopt.* ...
																							 TTnonopt)],'ylim',[0 500]);
	legend('O(\Delta t + \Delta x^2)','O(\Delta t^2 + \Delta x^2)',...
				 ['\theta = ',num2str(thetanonopt),', \Delta t_{n+1} = 0.25\Delta t_n, \Delta x_{n+1} = 0.5\Delta x_n'],...
				 ['\theta = ',num2str(thetaopt),', \Delta t_{n+1} = 0.5\Delta' ...
					' t_n, \Delta x_{n+1} = 0.5\Delta x_n'],1);
end;