function scores=MyFuncScores(fun,theta,varargin)
%Calculates the numerical gradient of the function fun. Specifically
%designed when fun is a log likelihood function
%INPUTS:        
% fun:      The function
% theta:    The function parameters;
% varargin: The extra argument, fun may have;

% OUTPUTS:
% scores:   The numerical scores of the function
k=length(theta);
nn = size(theta,2);
if nn>1
    theta = reshape(theta,[k*nn,1]);
end
k=length(theta);    
h=max(abs(theta*eps^(1/3)),1e-8);
h=diag(h);

[LLF,like]=feval(fun,theta,varargin{:});

t=length(like);

LLFp=zeros(k,1);
LLFm=zeros(k,1);
likep=zeros(t,k);
likem=zeros(t,k);
for i=1:k
    thetaph=theta+h(:,i);
    [LLFp(i),likep(:,i)]=feval(fun,thetaph,varargin{:});
    thetamh=theta-h(:,i);
    [LLFm(i),likem(:,i)]=feval(fun,thetamh,varargin{:});
end

scores=zeros(t,k);

h=diag(h);
for i=1:k
    scores(:,i)=(likep(:,i)-likem(:,i))./(2*h(i));
end