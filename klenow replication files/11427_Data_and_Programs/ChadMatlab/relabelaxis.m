function []=relabelaxis(vals,labs,whichaxis,axhandle);

% function []=relabelaxis(vals,labs,whichaxis);
%
% Assigns the labels in "labs" to the values in "vals" 
% to the axis label in "whichaxis"
%
% Note: labs can be a collection of cells or strmat or even just a string
%  e.g. labs={'$0' '$500k' '$1.0m' '$1.5m' '$2.0m' '$2.5m' '$3.0m'};

if ~exist('whichaxis')==1; whichaxis='y'; end;
if ~exist('axhandle'); axhandle=gca; end;
if size(labs,1)==1 & iscell(labs)~=1; labs=strmat(labs); end;
if whichaxis=='x';
    set(axhandle,'XTick',vals);
    set(axhandle,'XTickLabel',labs);
    % Also adjust XLim in case one of labels is outside the range:
    curlim=get(axhandle,'XLim');
    curlim(1)=min([curlim(1) min(vals)]');
    curlim(2)=max([curlim(2) max(vals)]');
    set(axhandle,'XLim',curlim);
else;
    set(axhandle,'YTick',vals);
    set(axhandle,'YTickLabel',labs);
    % Also adjust YLim in case one of labels is outside the range:
    curlim=get(axhandle,'YLim');
    %curlim(1)=min([curlim(1) min(vals)]');
    %curlim(2)=max([curlim(2) max(vals)]');
    curlim(1)=min(vals);
    curlim(2)=max(vals);
    set(axhandle,'YLim',curlim);
end;