function []=greygrid(yrs,vals);

% function []=greygrid(yrs,vals);
%
%  Creates grey gridlines on y axis at vals
definecolors;
hold on;
for i=1:length(vals);
    plot(yrs,vals(i)*ones(length(yrs)),'-','Color',mygrey,'LineWidth',2);
end;
set(gca,'Layer','top');  % Top put the grey line under the tick marks...
