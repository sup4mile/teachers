function [x,fval,flag]=fzerochad(f,x0,factor,NumTries,StopIfError,Verbose);

% Wrapper around fzero: x0=[xlow xhi].  Repeatedly update
%  by factor until find a sign change:  [xlow/factor(1) xhi*factor(2)].
%  Allows factor=[2 1.1] to update by different amounts in each direction.
%
%  f is the function to call.  We're looking for x s.t. f(x)=0

if exist('factor')~=1; factor=2; end;
if exist('NumTries')~=1; NumTries=5; end;
if exist('StopIfError')~=1; StopIfError=1; end; % Default is to stop...
if exist('Verbose')~=1; Verbose=0; end;
if length(factor)==1; factor=[factor factor]; end;
x00=x0;

%sign1=sign(f(x00(1)));
%sign2=sign(f(x00(2)));
f1=f(x00(1)); f2=f(x00(2)); 
i=1;
while sign(f1)==sign(f2) & i<NumTries;
  if x00(1)>0;  
      x00(1)=x00(1)/factor(1);
  else;
      x00(1)=x00(1)*factor(1);
  end;
  if x00(2)>0;
      x00(2)=x00(2)*factor(2);
  else;
      x00(2)=x00(2)/factor(2);
  end;

  f1=f(x00(1)); f2=f(x00(2)); 
  if Verbose;
      fprintf('  low=%10.5f  f(low)=%10.5f  | hi=%10.5f  f(hi)=%10.5f\n',[x00(1) f1 x00(2) f2]);
  end;
  i=i+1;
end;
%if sign1==sign2; disp 'No sign change found in fzerochad. Stopping...'; keyboard; 
if sign(f1)==sign(f2) | isinf(f1) | isinf(f2) | ~isreal(f1) | ~isreal(f2) | isnan(f1.*f2); %disp 'No sign change found in fzerochad. Assigning a NaN...';
  x=NaN; fval=NaN; flag=-1;
  if StopIfError; disp 'No sign change found in fzerochad or Inf/Imag problem. Stopping...'; keyboard; end;
else; 
  [x,fval,flag]=fzero(f,x00);
end;

