function gr=growthrate4(x,XYears,GrowthYears);

% growthrate.m
%
%  XYears = [1 1000 1500 1600 1700 1800 1900 2000]'
%
%  GrowthYears=[1 1000;
%       1000 1500;
%	1500 2000];
%
%  Will search XYears to find the appropriate indexes...


gr=zeros(size(GrowthYears,1),1);
for i=1:size(GrowthYears,1);
    indx1=find(XYears==GrowthYears(i,1));
    indx2=find(XYears==GrowthYears(i,2));
    T=GrowthYears(i,2)-GrowthYears(i,1);
  gr(i)=1/T*log(x(indx2)/x(indx1));
end;