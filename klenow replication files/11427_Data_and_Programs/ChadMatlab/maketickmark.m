function []=maketickmark(x,y,tsize,lwidth,fcolor);

% Draw a tick mark of length tsize (units are percent, e.g. .03) and linewidth lwidth centered
% at (x,y)

if exist('tsize')~=1; tsize=10; end;
if exist('lwidth')~=1; lwidth=1; end;
if exist('fcolor')~=1; fcolor='k'; end;
hold on;
ax=axis;
step=tsize/2*(ax(4)-ax(3));
plot([x x],[y-step y+step],'-','Color',fcolor,'LineWidth',lwidth);
