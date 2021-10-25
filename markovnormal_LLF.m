function [LLF, LL, theta_s0, theta_s1, pf1, pf2, pplus1, pplus2] = markovnormal_LLF(kappa,data,thetabar)
% function [LLF, LL, theta_s0, theta_s1, pf1, pf2, pplus1, pplus2] = markovnormal_LLF(kappa,data,thetabar)
% The negative copula log-likelihood of the Normal copula 
% with time-varying regime switching dependence
%
% Fri, 08 Mar, 2013.
%
% Osvaldo Candido Silva Filho
%
% INPUTS: kappa = [constant for state 0,
%                  constant for state 1,
%                  beta - autoregressive parameter,
%                  alpha - forcing variable parameter,
%                  probability p,
%                  probability q];
%         data = [U V]; uniform(0,1) margins 
%		  thetabar, results from estimation of unconditional model
%
% OUTPUTS: LLF negative log-likelihood value
%          LL scores
%          theta_s0,theta_s1 Parameter dynamic on each
%          state
%          pf1, pf2     Predicted probabilities
%          pplus1, pplus2  Filtered probabilities
%
% http://candidof.wordpress.com

n = size(data,1);

u = data(:,1);
v = data(:,2);

x = norminv(data(:,1),0,1);
y = norminv(data(:,2),0,1);

%Initials

w0=kappa(1);
w1=kappa(2);
b=kappa(3);
a=kappa(4);
p=(exp(kappa(5)))./(1+exp(kappa(5)));
q=(exp(kappa(6)))./(1+exp(kappa(6)));

theta_s0=ones(n,1);
theta_s1=ones(n,1);

theta_s0(1)=thetabar;
theta_s1(1)=thetabar;

p1t = .5*ones(n,1);
p2t = .5*ones(n,1);
p1 = .5*ones(n,1);
p2 = .5*ones(n,1);

%Initials for Pr[s0=j|I0], j=0,1.

p1t(1)=(1-p)./(2-p-q); %regime 1
p2t(1)=1-p1t(1);       %regime 0

for i=1:2;
for jj = 2:n;
   
   %regime0
   
    if jj<=10
        psi_s0 = w0 + a*mean(x(1:jj-1).*y(1:jj-1)) + b*theta_s0(jj-1);
    else
        psi_s0 = w0 + a*mean(x(jj-10:jj-1).*y(jj-10:jj-1)) + b*theta_s0(jj-1);
    end

    theta_s0(jj) = 1.998/(1+exp(-psi_s0))-0.999;

   %regime1
    if jj<=10
        psi_s1 = w1 + a*mean(x(1:jj-1).*y(1:jj-1)) + b*theta_s1(jj-1);
    else
        psi_s1 = w1 + a*mean(x(jj-10:jj-1).*y(jj-10:jj-1)) + b*theta_s1(jj-1);
    end

    theta_s1(jj) = 1.998/(1+exp(-psi_s1))-0.999;

 
%Filtered probabilities

    if jj==2
        p1(jj)=p.*p1t(jj-1)+(1-p).*(1-p1t(jj-1));
        p2(jj)=q.*(1-p1t(jj-1))+(1-q).*p1t(jj-1);
        
    elseif jj<n %Updating Pr[st=j|It], j=0,1.
        
        c1 = NormalCopula_pdf(u(jj-1),v(jj-1),theta_s1(jj-1));
        
        c1=c1.*p1(jj-1);
       
        c0 = NormalCopula_pdf(u(jj-1),v(jj-1),theta_s0(jj-1));
        
        c0=c0.*p2(jj-1);
      
        p1t(jj-1)=c1./(c0+c1);
        p2t(jj-1)=c0./(c0+c1);
        
        p1(jj)=p.*p1t(jj-1)+(1-p).*(1-p1t(jj-1));
        p2(jj)=q.*(1-p1t(jj-1))+(1-q).*p1t(jj-1);
    
    elseif jj==n && i==1
        
        pplus1=p1t;
        pplus2=p2t;
        
        sm1=pplus1;
        sm2=pplus2;
        
        for t=(n-2):-1:2
            sm1(t)=((p*pplus1(t)*sm1(t+1))/(p*pplus1(t)+(1-q)*pplus2(t)))+((1-p)*pplus1(t)*sm2(t+1))/((1-p)*pplus1(t)+q*pplus2(t));
            sm2(t)=1-sm1(t);
        end
   
        p1t(1)=sm1(1);
    
    end
end
end

pf1=p1;
pf2=p2;

pplus1=p1t;
pplus2=p2t;

     LL = -1*log((NormalCopula_pdf(u,v,theta_s1)).*pf1+(NormalCopula_pdf(u,v,theta_s0)).*pf2);
     LLF = sum(LL);

if isreal(kappa)==0;
   LLF = 1e6;
elseif isreal(LLF)==0;
   LLF = 1e7;
elseif isnan(LLF)
   LLF = 1e8;
elseif isinf(LLF)
   LLF = 1e9;
end