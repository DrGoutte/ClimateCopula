function [LLF, LL, TAU1_s0,TAU2_s0, TAU1_s1,TAU2_s1, pf1, pf2, pplus1, pplus2] = markovsjc_LLF(kappa,data,thetabar)
% function [LLF, LL, TAU1_s0,TAU2_s0, TAU1_s1,TAU2_s1, pf1, pf2, pplus1, pplus2] = markovsjc_LLF(kappa,data,thetabar)
% The negative copula log-likelihood of the SJC copula 
% with time-varying regime switching dependence
%
% Fri, 08 Mar, 2013.
%
% Osvaldo Candido Silva Filho
%
% INPUTS: kappa = [constant for tau1 on state 0,
%                  constant for tau2 on state 0,
%                  constant for tau1 on state 1,
%                  constant for tau2 on state 1,
%                  beta1 - autoregressive parameter for tau1 dynamics,
%                  alpha1 - forcing variable parameter for tau1 dynamics,
%                  beta2 - autoregressive parameter for tau2 dynamics,
%                  alpha2 - forcing variable parameter for tau2 dynamics,
%                  probability p,
%                  probability q,
%         data = [U V]; uniform(0,1) margins 
%		  thetabar = [tau1,tau2], results from estimation of unconditional model
%
% OUTPUTS: LLF negative log-likelihood value
%          LL scores
%          TAU1_s0,TAU2_s0, TAU1_s1,TAU2_s1  Parameter dynamic on each
%          state
%          pf1, pf2     Predicted probabilities
%          pplus1, pplus2  Filtered probabilities
%
% http://candidof.wordpress.com

n = size(data,1);
u = data(:,1);
v = data(:,2);

%Initials

w01=kappa(1);
w02=kappa(2);
w11=kappa(3);
w12=kappa(4);
b1=kappa(5);
b2=kappa(6);
a1=kappa(7);
a2=kappa(8);
p=(exp(kappa(9)))./(1+exp(kappa(9)));
q=(exp(kappa(10)))./(1+exp(kappa(10)));

TAU1_s0=ones(n,1);
TAU2_s0=ones(n,1);
TAU1_s1=ones(n,1);
TAU2_s1=ones(n,1);

TAU1_s0(1)=thetabar(1);
TAU2_s0(1)=thetabar(2);
TAU1_s1(1)=thetabar(1);
TAU2_s1(1)=thetabar(2);

p1t = .5*ones(n,1);
p2t = .5*ones(n,1);
p1 = .5*ones(n,1);
p2 = .5*ones(n,1);

%Initials for Pr[s0=j|I0], j=0,1.


p1t(1)=(1-p)./(2-p-q); %regime 1
p2t(1)=1-p1t(1);       %regime 0

for i = 1:2;
for jj = 2:n;
   
   %regime0
   if jj<=10
      psi1_s0 = w01 + b1*TAU1_s0(jj-1) + a1*(mean(abs(u(1:jj-1)-v(1:jj-1))));
      psi2_s0 = w02 + b2*TAU2_s0(jj-1) + a2*(mean(abs(u(1:jj-1)-v(1:jj-1))));
   else
      psi1_s0 = w01 + b1*TAU1_s0(jj-1) + a1*(mean(abs(u(jj-10:jj-1)-v(jj-10:jj-1))));
      psi2_s0 = w02 + b2*TAU2_s0(jj-1) + a2*(mean(abs(u(jj-10:jj-1)-v(jj-10:jj-1))));
   end
   
   TAU1_s0(jj) = 0.998/(1+exp(-psi1_s0)) + 0.001;		
   TAU2_s0(jj) = 0.998/(1+exp(-psi2_s0)) + 0.001;
   
   %regime1
   if jj<=10
      psi1_s1 = w11 + b1*TAU1_s1(jj-1) + a1*(mean(abs(u(1:jj-1)-v(1:jj-1))));
      psi2_s1 = w12 + b2*TAU2_s1(jj-1) + a2*(mean(abs(u(1:jj-1)-v(1:jj-1))));
   else
      psi1_s1 = w11 + b1*TAU1_s1(jj-1) + a1*(mean(abs(u(jj-10:jj-1)-v(jj-10:jj-1))));
      psi2_s1 = w12 + b2*TAU2_s1(jj-1) + a2*(mean(abs(u(jj-10:jj-1)-v(jj-10:jj-1))));
   end
   
   TAU1_s1(jj) = 0.998/(1+exp(-psi1_s1)) + 0.001;		
   TAU2_s1(jj) = 0.998/(1+exp(-psi2_s1)) + 0.001;
   
  
%Filtered probabilities
 
    if jj==2
        
        p1(jj)=p.*p1t(jj-1)+(1-p).*(1-p1t(jj-1));
        p2(jj)=q.*(1-p1t(jj-1))+(1-q).*p1t(jj-1);
        
    elseif jj<n %Updating Pr[st=j|It], j=0,1.
        
        c1=sym_jc_pdf(u(jj-1),v(jj-1),TAU1_s1(jj-1),TAU2_s1(jj-1));
        
        c1=c1.*p1(jj-1);
       
        c0=sym_jc_pdf(u(jj-1),v(jj-1),TAU1_s0(jj-1),TAU2_s0(jj-1));
            
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
        jj=2;  
   end
end
end

pf1=p1;
pf2=p2;

pplus1=p1t;
pplus2=p2t;

    LL1 = sym_jc_pdf(u,v,TAU1_s0,TAU2_s0).*pf2+sym_jc_pdf(u,v,TAU1_s1,TAU2_s1).*pf1;
    LL = -1.*log(LL1);
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