function [b,se,rsq]=addolsline(x,y,color,addtext,weights);

% For weights, we weight each observation by sqrt(weight)
% in computing the regression

  if exist('color')~=1; color='k'; end;
  if exist('addtext')~=1; addtext=0; end;
  if exist('weights')~=1; weights=ones(length(x),1); end;
  hold on;
  
  smpl=~isnan(x.*y.*weights);
  xx=x(smpl).*sqrt(weights(smpl));
  yy=y(smpl).*sqrt(weights(smpl));
  myones=ones(length(xx),1).*sqrt(weights(smpl));
  [b tstat sigma2 vcv rsq]=lstiny(yy,[myones xx]);
  se = b(2)/tstat(2);

  N=length(x(smpl));
  u     = y(smpl)-[ones(N,1) x(smpl)]*b;
  ybar  = 1/N*(ones(N,1)'*y(smpl));
  rsqUW   = 1-(u'*u)/(y(smpl)'*y(smpl) - N*ybar^2);

  % Note rsq is the RSquared from the weighted regression
  % while rsqUW is from the unweighted fit. 
  %
  %  Stata does not report either of these. I'm guessing that Stata
  %  reports the R2 from the weighted regression after clustering???
  
  bstr=sprintf('OLS Slope = %7.3f\n',b(2));
  sstr=sprintf(' Std. Err.    = %7.3f\n',se);
  rstr=sprintf('        R^2      = %7.2f\n',rsqUW);

  yfit=[ones(length(x),1) x]*b;
  [xmin,imin]=min(x);
  [xmax,imax]=max(x);
  plot([xmin xmax],yfit([imin imax]),'-','Color',color,'LineWidth',1);

  if addtext==1;
    putstr(bstr);
    putstr(sstr);
    putstr(rstr);
  end;
  