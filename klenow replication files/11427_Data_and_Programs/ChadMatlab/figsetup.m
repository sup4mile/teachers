% figsetup.m
%
%  A simple set of code to set the default properties of the figure
%   -- especially its size and position on the screen

clf;
set(gcf,'Outerposition',[321         194        1080         812]);
%set(gcf,'paperpositionmode','auto')
hold on;
printopt; % To be sure -depsc gets passed!