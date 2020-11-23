% PACKR.m    Deletes all rows with missing values from a matrix
%	     Just like Gauss's packr.

function x = packr(x);
[n k] = size(x);
if k==1;
	x(~isfinite(x))=[];
else;
	x(any(~isfinite(x)'),:) = [];
end;
