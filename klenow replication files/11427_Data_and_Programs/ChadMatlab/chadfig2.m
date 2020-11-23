function []=chadfig2(xlab,ylab,shrinkfactor,dowait,FS);

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
%
%  9/12/16 -- Includes a trap to catch latex expressions like '\eta' and NOT upper them.

if exist('FS')~=1; FS=10; end;
%chadfig(upper(xlab),upper(ylab),shrinkfactor,dowait);

newylab=upper(ylab);
istart=strfind(ylab,'\');
if ~isempty(istart);
    indxspace=strfind(ylab,' ');
    if ~iscell(ylab);
        if isempty(indxspace) | all(indxspace<istart);
            iend=length(ylab);
        else;
            spc=indxspace(indxspace>istart);
            iend=spc(1);
        end;
        newylab(istart:iend)=ylab(istart:iend);
    end;
end;

newxlab=upper(xlab);
istart=strfind(xlab,'\');
if ~isempty(istart);
    indxspace=strfind(xlab,' ');
    if isempty(indxspace) | all(indxspace<istart);
        iend=length(xlab);
    else;
        spc=indxspace(indxspace>istart);
        iend=spc(1);
    end;
    newxlab(istart:iend)=xlab(istart:iend);
end;


chadfig(newxlab,newylab,shrinkfactor,dowait);
set(get(gca,'YLabel'),'FontName','Helvetica','FontSize',FS)
set(get(gca,'XLabel'),'FontName','Helvetica','FontSize',FS)
