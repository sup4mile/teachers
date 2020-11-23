function []=shadebetweenlines(t,y1,y2,ShadeColor);

% function []=shadebetweenlines(t,y1,y2);   5/29/14
%
%  Fills the area between y1(t) and y2(t) with ShadeColor 
%  Key example: drawing confidence intervals for a time series
%
%  Idea is here: http://stackoverflow.com/questions/6245626/matlab-filling-in-the-area-between-two-sets-of-data-lines-in-one-figure

T=[t; flipud(t)];
Y=[y1; flipud(y2)];              % create y values for out and then back
fill(T,Y,ShadeColor);                  