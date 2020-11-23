% VarianceExamples2.m  -- Conditional mean / cutoff
%
%  Suppose the top p percent of people from a distribution work.
%  How does the E[x|x>x(p)] / x(p) productivity change as p increases?
%   -- For Pareto (no change?), Frechet (similar), logNormal?


diarychad('VarianceExamples2');

theta   = 2.1      % Frechet/Pareto parameter
sigmaLN = 0.35   % from Variance_IncomeWages.m
N       = 1000000 % Number of draws

p=[0.0001  % doctors, lawyers
   0.0002
   0.0010
   0.0100
   0.050];

cutoffs=round((1-p)*N); % The absolute numbers.


ITER=250
PercentChange=zeros(ITER,4);
PercentChangeAll=zeros(ITER,4);

for iter=1:ITER;


i=rand(N,1); % Percentiles from a uniform, to invert to get draws
i=sort(i);   % Ascending order

% Frechet
x=(-log(i)).^(-1/theta);
for n=1:length(p);
    xkeep{n}=x(cutoffs(n):N);
    ratioF(iter,n)=mean(xkeep{n})/x(cutoffs(n));
end;

% Pareto
x=(1-i).^(-1/theta);
for n=1:length(p);
    xkeep{n}=x(cutoffs(n):N);
    ratioP(iter,n)=mean(xkeep{n})/x(cutoffs(n));
end;

% Log Normal
x=logninv(i,0,sigmaLN);
for n=1:length(p);
    xkeep{n}=x(cutoffs(n):N);
    ratioLN(iter,n)=mean(xkeep{n})/x(cutoffs(n));
end;

% Exponential: p=F(x)=1-exp(-lambda*x)  % Variance is 1/lambda^2
%  1-p = exp(-lambda x) ==> x=-1/lambda*log(1-p)
x=-1/lambda*log(1-i);
for n=1:length(p);
    xkeep{n}=x(cutoffs(n):N);
    ratioE(iter,n)=mean(xkeep{n})/x(cutoffs(n));
end;

v=[ratioP(iter,:); ratioF(iter,:); ratioLN(iter,:); ratioE(iter,:)];
PercentChange(iter,:)=v(:,2)./v(:,1)-1;
PercentChangeAll(iter,:)=v(:,end)./v(:,1)-1;

end; %ITER

vmean=[mean(ratioP); mean(ratioF); mean(ratioLN); mean(ratioE)];

rtle={'Pareto','Frechet','LogNormal','Exponential'};
tle=sprintf('%7.4f',p);
tle(1)=[];
disp '---------------------------------------------------';
disp 'Ratio: E[x|x>x(p)]/x(p) of productivity for the Top p Percentiles';
fprintf('   Averaged across %5.0f simulations\n',ITER);
disp '---------------------------------------------------';
disp ' ';
cshow(rtle,vmean,'%10.4f',tle);
disp ' '; 
fprintf('Average Percent change in ratio: %6.4f to %6.4f\n',p(1:2));
cshow(rtle,mean(PercentChange)','%10.4f');
disp ' ';
fprintf('Average Percent change in ratio: %6.4f to %6.4f\n',[p(1) p(end)]);
cshow(rtle,mean(PercentChangeAll)','%10.4f');


diary off;
