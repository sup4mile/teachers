function h=plotlogline(x,y,sym,axlabel,color,LineWidth);

%  Plots y vs x with y on a log scale (natural log). Line plot.
%
%  sym is a symbols argument: '-' or '--' for example.
%  axlabel is the string to use to label the y axis.

definecolors;
if exist('sym')~=1; sym='-';end;
if exist('axlabel')~=1; axlabel=[];end;
if exist('color')~=1; color=myblue; elseif isempty(color); color=myblue; end;
if exist('LineWidth')~=1; LineWidth=LW; end;

h=[];
h=plot(x,log(y),sym,'Color',color,'LineWidth',LineWidth);

% Now fix the log scale for the y axis.
%val=[1000 2000 4000 8000 16000 32000]';
%labs=strmat('1000 2000 4000 8000 16000 32000');
if isempty(axlabel);
   curlabs=get(gca,'YTickLabel');
   %newnum=exp(str2num(curlabs));cellfun(@str2num,curlabs)
   newnum=exp(cellfun(@str2num,curlabs));
   newlabs=num2str(newnum,'%6.0f');
   if any(delta(str2num(newlabs))==0);
%   if size(newlabs,2)==1; 		% Add another significant digit
      newlabs=num2str(newnum,'%6.1f');
      set(gca,'YTick',log(1/10*round(newnum*10)));
   else;
      set(gca,'YTick',log(round(newnum)));
   end;
   set(gca,'YTickLabel',newlabs);
%  Old error:   set(gca,'YTick',str2num(curlabs));  %But there is rounding!!!!
else;
   if size(axlabel,1)==1; labs=strmat(axlabel); else; labs=axlabel; end;
   logval=log(str2num(axlabel));
   set(gca,'YTick',logval);
   set(gca,'YTickLabel',labs);
   % Also adjust YLim in case one of labels is outside the range:
   curlim=get(gca,'YLim');
   curlim(1)=min([curlim(1) min(logval)]');
   curlim(2)=max([curlim(2) max(logval)]');
   set(gca,'YLim',curlim);
end;
