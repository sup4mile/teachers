function []=gridlines(whichaxis,linewidth,transparency);

Not finished...

% function []=gridlines(whichaxis,transparency);
%
%  whichaxis = {'y', 'x', 'b' for both}
%  transparency = .5 to lighten grid lines
%
%  Currently only implemented for Y axis

if exist('linewidth')~=1; linewidth=1; end;
if exist('transparency')~=1; transparency=1; end;

if whichaxis~='x';
    yticks=get(gca,'YTicks');
    xlims=get(gca,'XLim');
    hold on;
    for i=1:length(yticks);
        plot(xlims,[yticks(i) yticks(i)],'-',
end;