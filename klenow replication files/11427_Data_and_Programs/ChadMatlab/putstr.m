% putstr.m  -- []=putstr(blah);
%
%  Uses mergtex(ntex(blah)) to place string.

function []=putstr(blah,color,fontsize,fontname);
if exist('color')~=1; color='k'; end;
if exist('fontname')~=1; fontname='Times'; end;
if exist('fontsize')~=1; fontsize=14; end;
disp(['Place the string: ' blah]);
%mergtex(ntex(blah));
gtext(blah,'FontName',fontname,'Color',color,'FontSize',fontsize);
