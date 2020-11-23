function []=chadfigyy(axhandle,newlabel,changey,chadfigversion,FS);

% Adds a label to the second axis in a plotyy plot. Call after chadfig...
% Pass axhandle ==> ax=plotyy  axhandle=ax(2)

y2handle=get(axhandle,'YLabel');
if exist('changey')~=1; changey=1.7; end;
if exist('chadfigversion')~=1; chadfigversion=2; end;
if exist('FS')~=1; FS=10; end;

pos=get(y2handle,'Position');
curaxis=axis(axhandle);
step=curaxis(4)-curaxis(3);
newaxis=curaxis;
newloc=curaxis(4)+.04*step*changey;  % .03 instead of 0.1
                                     
%axis(newaxis);
pos(2)=newloc;
set(y2handle,'Position',pos);
set(y2handle,'HorizontalAlignment','right');
if chadfigversion==2;
    set(y2handle,'String',upper(newlabel),'Rotation',0,'FontName','Helvetica','FontSize',FS,'FontWeight','bold');
else;
    set(y2handle,'String',newlabel,'Rotation',0,'FontWeight','bold');
end;