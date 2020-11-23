function means=weightedmean(X,weights,ReturnMedian);

% weightedmean.m  11/19/99
%
%                  function means=weightedmean(X,weights,ReturnMedian);
%
%  Function to compute the weighted mean of the columns of a matrix. 
%     X         = The matrix
%   weights = The weights (need not sum to 1; we'll normalize).
%
%   If ReturnMedian==1 then return the median of unweighted observations instead.

w=weights/sum(weights);
Xw=mult(X,w); 				% Multiply each column by w
means=sum(Xw); 				% And then add up...

if exist('ReturnMedian');   % Used in my Compustat work in IdeaPF
    if ReturnMedian==1;
        means=median(X);
    end;
end;