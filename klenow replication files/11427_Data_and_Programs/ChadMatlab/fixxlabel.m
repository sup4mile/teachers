function []=fixxlabel(relative_offset);

% Fixing the spacing between xlabel and xaxis
% See here: http://stackoverflow.com/questions/14966770/distance-between-axis-label-and-axis-in-matlab-figure
if exist('relative_offset')~=1;
    relative_offset = 2;
end;
xh = get(gca,'XLabel'); % Handle of the x label
set(xh, 'Units', 'Normalized')
pos = get(xh, 'Position');
set(xh, 'Position',pos.*[1,relative_offset,1])
