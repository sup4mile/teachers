function x=generatePareto(N,eta,x0);

% function x=generatePareto(N,eta,x0);
%  5/30/14
%
%   Generate N draws from a pareto distribution with minimum value x0 and
%   pareto inequality parameter eta
%
%   F(x) = 1 - (x/x0)^(-1/eta)
%
%   Idea: F(.) is uniform and invert:
%    
%     x = x0*(1-F)^(-eta)

F=rand(N,1);
x=x0*(1-F).^(-eta);

