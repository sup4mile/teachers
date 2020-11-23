function []=shade(start,finish,colorstr,yaxisadj);

% function []=shade(start,finish,colorstr);
%
%  start and finish are Nx1 vectors of starting and ending years.
%  The function shades between the start and finish pairs using colorstr

if ~exist('colorstr'); colorstr='y'; end;  % default is yellow
if ~exist('yaxisadj'); yaxisadj=0; end;
curax=axis;
y0=curax(3)+yaxisadj*(curax(4)-curax(3)); % Added 9/17/15 to start just above the axis
y=[y0 curax(4) curax(4) y0];
hold on;
for i=1:length(start);
  x=[start(i) start(i) finish(i) finish(i)];
  fill(x,y,colorstr);
end;
  
% Now, prevent the shading from covering up the lines in the plot.  
%h = findobj(gca,'Type','line');
%set(h,'EraseMode','xor');
% EraseMode no longer works in Matlab 2015...
kids = get(gca,'Children');        %# Get the child object handles
set(gca,'Children',flipud(kids));  %# Set them to the reverse order

h = findobj(gca,'Type','patch');
set(h,'EdgeColor','none');

% This last one makes the tick marks visible
if yaxisadj==0;
    set(gca, 'Layer', 'top')
end;
