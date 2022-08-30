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
% FINAL PART(computation of the solution depending on the RHS fi)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
       function fi=NSE_ADI_step(ami,api,alph,xs2,fi)

       [m,n]=size(fi);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %
% Resolution of the tridiagonal system                             %
%  ami(i)*X(j,i-1)+aci(i)*X(j,i)+api(i)*X(i+1)=fi(j,i), i=1,n      %
%   solution stored  in fi                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%***------------------------------------------------------------
%***   Forward substitution
%***------------------------------------------------------------

      for i=2:n
	fi(:,i)= fi(:,i)-ami(i)*fi(:,i-1);
      end


%***------------------------------------------------------------
%***  Backward substitution  -->	dudx1(i)=X1(i)
%                            =gamma(n), i=n
%                            =gamma(i)-[c(i)/beta(i)]*X1(i+1)
%***------------------------------------------------------------

        fi(:,n)=fi(:,n)*alph(n);
%

        for i=n-1:-1:1
            fi(:,i)=fi(:,i)*alph(i)-api(i)*fi(:,i+1);
        end
%***------------------------------------------------------------
%***   Final solution    	-->  v2_x = X*
%                              dudx1 = X = [X1] - [X*] [X2]
%***------------------------------------------------------------
	  xs1=zeros(m,1);
          xs1=(fi(:,1)+fi(:,n))/(1.d0+xs2(1)+xs2(n));

%
        for  i=1:n
          fi(:,i)=fi(:,i)-xs1*xs2(i);
        end