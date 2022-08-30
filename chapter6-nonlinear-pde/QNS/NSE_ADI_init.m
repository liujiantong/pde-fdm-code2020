%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    An Introduction to Scientific Computing          %%%%%%%
%%%%%%%    I. Danaila, P. Joly, S. M. Kaber & M. Postel     %%%%%%%
%%%%%%%                 Springer, 2005                      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                %
%  Solves simultaneously m linear systems                        %
% ami(i)*X(j,i-1)+aci(i)*X(j,i)+api(i)*X(j,i+1)=fi(j,i),         %
%               i=1:n                                            %
%  for    each  j=1:m                                            %
%  periodicity condition  X(j,1)=X(j,n)                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INITIALIZATION PART (computations not depending on the RHS)   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
       function [ami,api,alph,xs2]=NSE_ADI_init(ami,aci,api)

       n=length(ami);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modification of the matrix (1st and last coeff. on the diag.)  %
%         aci(1)=aci(1)-ami(1)                                   %
%         aci(n)=aci(n)-api(n)                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
         aci(1)=aci(1)-ami(1);
         aci(n)=aci(n)-api(n);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of the RHS for the calculation of  X2           %
%         xs2(1)=ami(1)                                          %
%         xs2(i) =0                                              %
%         xs2(n)=api(n)                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
	       xs2=zeros(1,n);
         xs2(1)=ami(1);    
         xs2(n)=api(n);    

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %
% Resolution of the tridiagonal system                             %
% ami(i)*X(i-1)+aci(i)*X(i)+api(i)*X(i+1)=xs2(i),                  %
%      i=1,n                                                       %
%   the solution is stored in  xs2                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%--> Thomas algorithm 
%-->                 alph(i)=1/beta(i)
%                     ami(i)=ami(i)/beta(i-1)
%                     api(i)=api(i)/beta(i)
          alph=zeros(1,n);
          alph(1)=1.d0/aci(1);
        for i=2:n
          alph(i)=1.d0/(aci(i)-ami(i)*api(i-1)*alph(i-1));
          ami(i)=ami(i)*alph(i-1);
          api(i-1)=api(i-1)*alph(i-1);
	end

%
%-->resolution of the system [M*] X2 = v1_x
%   attention,  aci(i)=beta(i)*gamma(i)

          aci(1)=xs2(1);
        for i=2:n
	  aci(i)=xs2(i)-ami(i)*aci(i-1);
        end

%                       xs2(i) = X2(i) 
	xs2(n)=aci(n)*alph(n);
      for i=n-1:-1:1
        xs2(i)=aci(i)*alph(i)-api(i)*xs2(i+1);
      end
