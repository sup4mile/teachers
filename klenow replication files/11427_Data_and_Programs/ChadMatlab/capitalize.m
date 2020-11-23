function newstr=capitalize(str);

% INTEL COMPANY -> Intel Company
% http://www.mathworks.com/matlabcentral/answers/107307-function-to-capitalize-first-letter-in-each-word-in-string-but-forces-all-other-letters-to-be-lowerc

YesCell=0;
if iscell(str);
  YesCell=1;
    
  % Create strmat version of text for cshow
  a=str(1);
  for i=2:length(str);
    a=[a '#' str(i)];
  end;
  str=strmat(cell2mat(a),'#');
end;

N=size(str,1)
for i=1:N;
    stri=lower(str(i,:));
    idx=regexp([' ' stri],'(?<=\s+)\S','start')-1;
    stri(idx)=upper(stri(idx));
    newstr(i,:)=stri;
end;

if YesCell;
    newstr=cellstr(newstr);
end;