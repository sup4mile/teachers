function []=testfzero(fname,xvals);

% []=testfzero(fname,xvals);     Example:  testfzero(fun,(.1:.01:.5))
%    fname = function we are trying to find a zero of
%    xvals = e.g. (-.5:.01:.5) -- range of x's we are considering.
%
%  Plots a graph of fname(xvals) 

for i=1:length(xvals);
    fval(i)=feval(fname,xvals(i));
end;
clf;
plot(xvals,fval);