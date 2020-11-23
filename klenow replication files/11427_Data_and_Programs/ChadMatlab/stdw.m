
function s=stdw(x,w,flag);

% stdw.m Computed weighted standard deviation. Just calls std(x,w,flag)
%   except that std(.) just wants a vector of weights. This stdw function
%   allows the weights to differ by column, i.e. w is a matrix of weights
%   with the same dimension as x.

if exist('flag')~=1; flag='includenan'; end;

if size(w,2)==1;
    s=std(x,w,flag);
else;
    for i=1:size(w,2);
        s(i)=std(x(:,i),w(:,i),flag);
    end;
end;