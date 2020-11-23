function [avg,stdev]=weightedaverage(x,weight,keepit);

% weightedaverage.m   [avg,stdev]=weightedaverage(x,weight,keepit);
%
% x = NxK matrix to average (columnwise)
% weight = weighting variable, Nx1, e.g. population - need not sum to one (that's what keepit is for)
% keepit = logical Nx1 with "1" corresponding to observations to keep.
%
% As a bonus, also computes the stdev of the weighted data. See formula here
% (only implemented for a vector at this point)
%
% http://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf

if exist('keepit')~=1; keepit=ones(size(x,1),1); end;
totweight=sum(weight(keepit));

if size(x,2)~=1;
    stdev=[];
    avg=sum(mult(x(keepit,:),weight(keepit)/totweight));
else;
    avg=sum(x(keepit).*weight(keepit)/totweight);
    N=length(x(keepit));
    s2=sum(weight(keepit).*(x(keepit)-avg).^2);
    stdev=sqrt(N*s2 / ((N-1)*totweight));
end;
