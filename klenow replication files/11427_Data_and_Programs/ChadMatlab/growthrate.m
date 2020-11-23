function gr=growthrate(x,yrs);

% growthrate.m
%
%  yrs=[1 10;
%       1 5;
%	6 10];
%
%  x=TxN vector  ==> compute the average growth rate of x over the yrs periods.

gr=zeros(size(yrs,1),size(x,2));
for i=1:size(yrs,1);
    if size(x,2)>1;
        gr(i,:)=mean(delta(log(x(yrs(i,1):yrs(i,2),:))));
    else;
        gr(i)=log(x(yrs(i,2))/x(yrs(i,1)))/(yrs(i,2)-yrs(i,1));
    end;
end;