% LSweight.m 
%
%  b=inv(x'wx)*x'wy  where w is diagonal matrix of weights

function [beta,tstat,sigma2,vcv,rsq] = lsweight(y,x,w);

[N K] = size(x);
xxinv = inv(x'*w*x);
beta  = xxinv*x'*w*y;
u     = y-x*beta;
sigma2=u'*u/(N-K);
stdest=sqrt(sigma2);
vcv   =sigma2*xxinv;
se    =sqrt(diag(vcv));
tstat = beta./se;

????
ybar  = 1/N*(ones(N,1)'*y);
rsq   = 1-(u'*w*u)/(y'*w*y - N*ybar^2);


See addolsline.m
