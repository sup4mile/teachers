function ydata=fixmissing(x,missval,NotStartEnd);

% function xnew=fixmissing(x,missval,NotStartEnd);
%
%   x = matrix to be "fixed" -- will fix each column
%   missval = boolean of the entries to be fixed (optional
%     if missval is not included, use isnan(x).
%
%   xnew = replace with first not missing / last not missing / average of two surrounding
%
%    If "NotStartEnd" is absent or zero, missing values at beginning or end
%    will be replaced with first/last value.
%  
%    If NotStartEnd==1, then leaves missing values alone at beginning and end
%    If NotStartEnd==2, then leaves missing values alone at beginning only
%    If NotStartEnd==3, then leaves missing values alone at end only

% if exist('NotStartEnd')==0;
%   NotStartEnd=0;
% end;

% if exist('missval')~=1; missval=isnan(x); end;
% indx=find(~missval);
% xnew=x;
% for i=1:length(x);
%     if i==1 & missval(i); xnew(i)=x(indx(1)); % first non-missing val
%     elseif i>indx(end) & missval(i); xnew(i)=x(indx(end)); % last non-missing val
%     elseif missval(i); xnew(i)
        
        
    
%     end;        
% end;


% function ydata=interplin3(x,missval,NotStartEnd)

% % interplin3
% %
% %   Simple linear interpolation -- replaces NaNs (or missval)
% %
% %    x=initial vector/matrix (will interpolate each column)
% %    missval = (optional) denotes "missing value", default is NaN
% %    Column by column
% %
% %    If "NotStartEnd" is absent or zero, missing values at beginning or end
% %    will be replaced with first/last value.
% %  
% %    If NotStartEnd==1, then leaves missing values alone at beginning and end
% %    If NotStartEnd==2, then leaves missing values alone at beginning only
% %    If NotStartEnd==3, then leaves missing values alone at end only

if exist('NotStartEnd')==0;
  NotStartEnd=0;
end;
if exist('missval')~=1; missval=isnan(x); end;

[T K]=size(x);
if T==1; x=x'; T=length(x); end;  % row vector ==> column vector

xdata=x;

for c=1:K;  % column by column
  x=xdata(:,c);
  if ~all(missval); % Skip column if all values are missing
    y=x;
    i0=find(~missval(:,c));
    iT=i0(length(i0)); 			% The last valid position
    i0=i0(1); 				% The first valid entry

    for i=i0:(iT-1);     % position in x
      if missval(i,c);
        notnan=find(~missval((i+1):iT,c));  % Next non NaN
        notnan=notnan(1)+i;
        y(i:(notnan-1))=(x(notnan)+x(i-1))/2;
        missval(i:(notnan-1),c)=0;
        %        keyboard
      end;
    end;

    ydata(:,c)=y;
    if NotStartEnd==0;
      ydata(1:(i0-1),c)=y(i0);  % Missing at beginning
      ydata((iT+1):end,c)=y(iT);  % Missing at end;
    elseif NotStartEnd==2;
        ydata((iT+1):end,c)=y(iT);  % Missing at end fixed
    elseif NotStartEnd==3;
        ydata(1:(i0-1),c)=y(i0);  % Missing at beginning fixed
    end;
  end;
end;

