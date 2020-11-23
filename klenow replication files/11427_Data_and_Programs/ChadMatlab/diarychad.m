function []=diarychad(fname,casename);

% function []=diarychad(fname,casename);
% 
% Deletes the file 'fname.casename.log' if it exits and opens it for diary.
% If fname.m exists, shows the help fname.m
%   (casename is optional)

if exist('casename');
    dname=[fname '_' casename];
else;
    dname=fname;
end;
if exist([dname '.log']); delete([dname '.log']); end;
diary([dname '.log']);
fprintf([dname '                 ' date]);
disp ' ';
disp ' ';
if exist([fname '.m']);
    eval(['help ' fname]);
end;
