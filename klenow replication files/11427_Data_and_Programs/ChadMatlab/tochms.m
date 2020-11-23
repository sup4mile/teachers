function t=tochms;

%t=datevec(toc./(60*60*24))
t=datestr(toc/(60*60*24), 'HH:MM:SS.FFF');
disp(['Timer checked after H:M:S = ' t]);