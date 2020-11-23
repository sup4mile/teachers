function newy=interplin4(y,x,newx)

% interplin4
%     function newy=interplin(y,x,newx)
%
%   Simple linear interpolation, assuming that x is monotone
%
%     y = Variable to be interpolated  -- Nx1
%     x = Values we know, corresponding to existing y -- Nx1
%   newx = x values for which we need y values (can include x)
%                    1960 10
%                    1965 15
%                    1967 14
%
%    newy=interplin(y,x,[1962; 1966]) will return [12; 14.5]
%
%   If newx is outside the range of x, return the endpoint
%   e.g. interplin(y,x,1968) returns 14.

newy=zeros(length(newx),1);
for i=1:length(newx);
   lessthan=find(x<=newx(i));
   greatthan=find(x>=newx(i));
   jlow=max(lessthan);
   jhigh=min(greatthan);
   if jlow==jhigh; 
      newy(i)=y(jlow); 
   elseif isempty(greatthan);
       newy(i)=y(jlow);
   elseif isempty(lessthan);
       newy(i)=y(jhigh);
   else;
      dx=x(jhigh)-x(jlow);
      dy=y(jhigh)-y(jlow);
      newy(i)=y(jlow)+dy/dx*(newx(i)-x(jlow));
   end;
end;

