function [CL, rhohat,nuhat] = bivt_tvp1_CL(theta,Zdata,rhobar)

T = size(Zdata,1);


x = tinv(Zdata(:,1),10);
y = tinv(Zdata(:,2),10);


rho = -999.99*ones(T,1);
nu = -999.99*ones(T,1);
rho(1) = rhobar(1);			
nu(1) = rhobar(2);
psi(1,:)=zeros(1,2);

for jj = 2:T
    if jj<=10
        psi1 = theta(1) + theta(2)*mean(x(1:jj-1).*y(1:jj-1)) + theta(3)*rho(jj-1);
        psi2 = theta(4) + theta(5)*mean(x(1:jj-1).*y(1:jj-1)) + theta(6)*nu(jj-1);
    else
        psi1 = theta(1) + theta(2)*mean(x(jj-10:jj-1).*y(jj-10:jj-1)) + theta(3)*rho(jj-1);
        psi2 = theta(4) + theta(5)*mean(x(jj-10:jj-1).*y(jj-10:jj-1)) + theta(6)*nu(jj-1);
    end
    
    psi(jj,:) = [psi1,psi2];
    rho(jj) = 1.998/(1+exp(-psi1))-0.999;		
    nu(jj) = (exp(psi2)/(1+exp(psi2)))*98 + 2;
end

rhohat = rho;  
nuhat = nu;


%for i=1:T

 %   rhom = [1 rho(i); rho(i) 1];

  %  CL(i) = copulapdf('t',[Zdata(:,1) Zdata(:,2)],rhom,nu(i));
%end
CL = tcopula_pdf(Zdata(:,1),Zdata(:,2),rho,nu);
CL = log(CL);

%for i=2:T
   
    
  %CL(i) = gammaln((nu(i)+2)./2) + gammaln(nu(i)./2) - 2*gammaln((nu(i)+1)./2) - 0.5*log(1-rho(i)^2);
  %CL(i) = CL(i) - (nu(i)+2)./2*log(1+(x(i).^2 + y(i).^2 - 2*rho(i)*x(i).*y(i))./(nu(i)*(1-rho(i)^2)));
  %CL(i) = CL(i) + (nu(i)+1)./2*log(1+x(i).^2/nu(i)) + (nu(i)+1)./2*log(1+y(i).^2/nu(i));
    
%end    
    
CL = sum(CL);
CL = -CL;