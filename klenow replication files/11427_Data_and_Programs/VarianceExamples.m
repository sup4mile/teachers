% VarianceExamples.m
%
%  Suppose the top p percent of people from a distribution work.
%  How does the variance of productivity change as p increases?
%   -- For Pareto (no change?), Frechet (similar), logNormal?


diarychad('VarianceExamples');

theta   = 2.1            % Frechet/Pareto parameter
sigmaLN = sqrt(0.35)   % from Variance_IncomeWages.m
lambda=1/sigmaLN       % Stdev(exp) = 1/lambda
N       = 1000000 % Number of draws

p=[0.0001  % doctors, lawyers
   0.0002
   0.0010
   0.0100
   0.050];

cutoffs=round((1-p)*N); % The absolute numbers.

ITER=250
VarianceRatio=zeros(ITER,4);
VarianceRatioAll=zeros(ITER,4);
VarianceRatio_log=zeros(ITER,4);
VarianceRatioAll_log=zeros(ITER,4);

for iter=1:ITER;



i=rand(N,1); % Percentiles from a uniform, to invert to get draws
i=sort(i);

% Frechet
x=(-log(i)).^(-1/theta);
for n=1:length(p);
    xkeep{n}=x(cutoffs(n):N);
    varFlog(iter,n)=variance(log(xkeep{n}));
    varF(iter,n)=variance((xkeep{n}));
end;

% Pareto
x=(1-i).^(-1/theta);
for n=1:length(p);
    xkeep{n}=x(cutoffs(n):N);
    varPlog(iter,n)=variance(log(xkeep{n}));
    varP(iter,n)=variance((xkeep{n}));
end;

% Log Normal
x=logninv(i,0,sigmaLN);
for n=1:length(p);
    xkeep{n}=x(cutoffs(n):N);
    varLNlog(iter,n)=variance(log(xkeep{n}));
    varLN(iter,n)=variance((xkeep{n}));
end;

% Exponential: p=F(x)=1-exp(-lambda*x)  % Variance is 1/lambda^2
%  1-p = exp(-lambda x) ==> x=-1/lambda*log(1-p)
x=-1/lambda*log(1-i);
for n=1:length(p);
    xkeep{n}=x(cutoffs(n):N);
    varElog(iter,n)=variance(log(xkeep{n}));
    varE(iter,n)=variance((xkeep{n}));
end;

varlog=[varPlog(iter,:); varFlog(iter,:); varLNlog(iter,:); varElog(iter,:)];
VarianceRatio_log(iter,:)=varlog(:,2)./varlog(:,1);
VarianceRatioAll_log(iter,:)=varlog(:,end)./varlog(:,1);

v=[varP(iter,:); varF(iter,:); varLN(iter,:); varE(iter,:)];
VarianceRatio(iter,:)=v(:,2)./v(:,1);
VarianceRatioAll(iter,:)=v(:,2)./v(:,1);

end; %ITER

vmeanlog=[mean(varPlog); mean(varFlog); mean(varLNlog); mean(varElog)];
vmean=[mean(varP); mean(varF); mean(varLN); mean(varE)];



rtle={'Pareto','Frechet','LogNormal','Exponential'};
tle=sprintf('%7.4f',p);
tle(1)=[];
disp '---------------------------------------------------';
disp 'Variances of productivity for the Top p Percentiles';
fprintf('   Averaged across %5.0f simulations\n',ITER);
disp '---------------------------------------------------';
disp ' ';
disp 'Variance of log x';
cshow(rtle,vmeanlog,'%10.4f',tle);
disp ' '; 
fprintf('Average ratio of variance: %6.4f to %6.4f\n',p(1:2));
cshow(rtle,mean(VarianceRatio_log)','%10.4f');
disp ' ';
fprintf('Average ratio of variance: %6.4f to %6.4f\n',[p(1) p(end)]);
cshow(rtle,mean(VarianceRatioAll_log)','%10.4f');

disp ' '; disp ' ';
disp 'Variance of x';
cshow(rtle,vmean,'%10.4f',tle);
disp ' '; 
fprintf('Average ratio of variance: %6.4f to %6.4f\n',p(1:2));
cshow(rtle,mean(VarianceRatio)','%10.4f');
disp ' ';
fprintf('Average ratio of variance: %6.4f to %6.4f\n',[p(1) p(end)]);
cshow(rtle,mean(VarianceRatioAll)','%10.4f');


diary off;
