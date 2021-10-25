function [CL, kappa] = Clayton_tvp1_CL(theta,data,kappabar)
%function [CL, kappa] = Galambos_tvp1_CL(theta,data,kappabar)
% 
% The negative copula log-likelihood of Galambos's copula with time-varying parameter 
%
%
% INPUTS: theta ;
%				data = [U V], a Tx2 matrix of data
%				thetabar = parameter of the Galambos copula without time-variation
%				exog, a Txk matrix of exogenous conditional correlation regressors

T = size(data,1);
u = data(:,1);
v = data(:,2);

kappa = -999.99*ones(T,1);
kappa(1) = kappabar;			
for jj = 2:T;
   if jj<=10
      psi1 = theta(1)+ theta(2)*kappa(jj-1) + theta(3)*(mean(abs(u(1:jj-1)-v(1:jj-1))));
   else
      psi1 = theta(1)+ theta(2)*kappa(jj-1) + theta(3)*(mean(abs(u(jj-10:jj-1)-v(jj-10:jj-1))));
   end
   kappa(jj) = 0.0001+psi1.^2;
end

CL = clayton_pdf(u,v,kappa);
CL = log(CL);
CL = sum(CL);
CL = -CL;

if isreal(theta)==0;
   CL = 1e6;
elseif isreal(CL)==0;
   CL = 1e7;
elseif isnan(CL)
   CL = 1e8;
elseif isinf(CL)
   CL = 1e9;
end