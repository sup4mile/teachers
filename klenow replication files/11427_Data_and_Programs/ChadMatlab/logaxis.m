function []=logaxis(

   if size(axlabel,1)==1; labs=strmat(axlabel); else; labs=axlabel; end;
   logval=log(str2num(axlabel));
   set(gca,'XTick',logval);
   set(gca,'XTickLabel',labs);
   % Also adjust XLim in case one of labels is outside the range:
   curlim=get(gca,'XLim');
   curlim(1)=min([curlim(1) min(logval)]');
   curlim(2)=max([curlim(2) max(logval)]');
   set(gca,'XLim',curlim);