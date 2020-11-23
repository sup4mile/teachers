function []=chadtimer;

% Assuming 'tic' has already been set, this function calls 'toc' and formats the output
% in Hrs/Mins/Secs

fprintf('Timer check: ');
t_hms = datevec(toc./(60*60*24));          % Express ‘t’ As Fraction Of A Day, Then Use ‘datevec’ To Convert
fprintf('%2.0f hours, %2.0f minutes, %2.0f seconds\n',t_hms(4:6));