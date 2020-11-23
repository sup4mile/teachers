function []=chadgraph2(xlab,ylab,shrinkfactor,dowait);

% chadfig2.m   3/12/04
%
%   function []=chadfig(xlab,ylab,shrinkfactor,dowait,changey);
%
%   changey=0 to leave the yaxis size unchanged.
%
%  Try to produce very high quality figures that look like those in my
%  growth book.
%
%   xlab=xlabel,  ylab=ylabel.  
%
%  Good ideas for nice figures, to implement before calling:
%
%   xlab='Productivity, $\bar{A}$' now works (latex interpreter added)
%   Actually, no, it doesn't work.  The problem is the latex interpreter gets
%   applied to the entire text, overriding bold, etc. and making it look weird.
%
%   xlab=upper(xlab)     %  Upper Case labels
%   set(gca,'LineWidth',1);
%   plot figure with LineWidth=1 as well.

chadgraph(upper(xlab),upper(ylab),shrinkfactor,dowait);
set(get(gca,'YLabel'),'FontName','Helvetica','FontSize',10)
set(get(gca,'XLabel'),'FontName','Helvetica','FontSize',10)
