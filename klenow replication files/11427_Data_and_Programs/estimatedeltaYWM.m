
% estimatedeltaYWM.m: Estimate delta from YWM using Education to identify phi
%   Passing "A" to use as instruments. Initially NaN, however...
%
%  1/28/19


function [dlta,Phi,S]=estimatedeltaYWM(eta,beta,theta,gam,sigma,Tbar,A);

%if exist('ShowData')~=1; ShowData=0; end;

global CaseName;
load(['CohortData_' CaseName '.mat']);

if exist('tauw')~=1; tauw=zeros(size(p)); end;
if exist('tauh')~=1; tauh=zeros(size(p)); end;

HOME=1;
Mkt=2:Noccs; % Market Occupation numbers

% Requires Ty/Tbar in order to use the wagebar equation.
Ty=ones(Noccs,1); % Implicitly assumed below, not passed to the functions.
Tbar=squeeze(Tbar(:,1,:)); % Pull WM ==> Noccs x Nyears

% Store results in Noccs x Decades matrices
W=zeros(Noccs,Nyears)*NaN;
PHI=W; S=W;
MWM=zeros(Nyears,1)*NaN;   % Nyears x 1
WageBarHomeY=zeros(Ngroups,Nyears)*NaN;
mu=1/theta*1/(1-eta);
etabar=eta^(eta/(1-eta));

for t=1:Nyears;

    disp ' '; disp ' ';
    disp '-------------------------------------------------------------';
    tlestr=sprintf('Estimating delta from YWM: %6.0f\n',CohortConcordance(t,1));
    disp(tlestr);
    disp '-------------------------------------------------------------';

    YoungCohort=CohortConcordance(t,2);
    wagebar=squeeze(WageBar(:,WM,YoungCohort,t)); 
    pt=squeeze(p(:,:,YoungCohort,t));
    pWM=pt(:,WM);

    % Use Education to get phi
    s=EducationYWM(:,t)/25; % 25 years is denominator
    phi=(1-eta)/beta*s./(1-s); % Noccs x 1
                               %phi(HOME)=sum(EarningsWeights_igt(Mkt,1,t).*phi(Mkt)); % Home -- use average of Market

    % Now adjust wagebar for (1-s) and Tbar/Ty (but these are assumed common across occs anyway)
    wagebarhat=wagebar.*(1-s).^(1/beta).*Tbar(:,t)./Ty;
    

    fmt='%8.3f %8.3f %10.1f %8.3f %8.3f %8.3f %8.1f %8.1f %8.3f';
    tle='s phi wagehat pWM';
    cshow(OccupationNames,[s phi wagebarhat pWM],fmt,tle);
    
    sfigure(1); figsetup;
    plotnamesym2(log(pWM),log(wagebarhat),ShortNames,9,[],.02,.00);
    chadfig2('p(i,t), log scale','WageBar (adj. for schooling, log)',1,0);
    makefigwide;
    title(num2str(Decades(t)));
    xt=[1/1024 1/256 1/64 1/16 1/4 1];
    set(gca,'XTick',log(xt));
    set(gca,'XTickLabel',strmat('1/1024# 1/256#  1/64#  1/16#   1/4#    1','#'));
    if t==1;
        print('-dpsc',['estimatedeltaYWM_' CaseName]);
    else;
        print('-dpsc','-append',['estimatedeltaYWM_' CaseName]);
    end;        
    [bcoeff,se,tstat,robse,trob]=ols2(log(wagebarhat),[ones(Noccs,1) log(pWM)],tlestr,'logwagebarhatYWM',strmat('Constant log(pYWM)'));
    disp 'Note: the slope is an estimate of delta / (theta(1-eta))';
    DeltaOLS(t,:)=[bcoeff(2)/mu robse(2)/mu trob(2)];
    
    W=[ones(Noccs,1) log(A(:,t))]; % Instruments
    instv=strmat('Constant log(A)');
    [bcoeff,se,tstat,robse,trob]=iv2(log(wagebarhat),[ones(Noccs,1) log(pWM)],W,tlestr,'logwagebarhatYWM',strmat('Constant log(pYWM)'),instv,1);
    DeltaIV(t,:)=[bcoeff(2)/mu robse(2)/mu trob(2)];
    
    Phi(:,t)=phi;
    S(:,t)=s;
end; % Looping over decades

disp ' ';
disp 'LONG DIFFERENCE: 1960 - 2010';

wagebar=[squeeze(WageBar(:,WM,7-1,1)) squeeze(WageBar(:,WM,7-6,6))];
s=[S(:,1) S(:,6)];
TTbar=[Tbar(:,1) Tbar(:,6)];
wagebarhat=wagebar.*(1-s).^(1/beta).*TTbar; % Since Ty=ones
dlogwagebarhat = log(wagebarhat(:,2)./wagebarhat(:,1));

pt=[squeeze(p(:,WM,7-1,1)) squeeze(p(:,WM,7-6,6))];
dlogp=log(pt(:,2)./pt(:,1));

sfigure(1); figsetup;
plotnamesym2(dlogp,dlogwagebarhat,ShortNames,9,[],.2,.05);
%chadfig2('dlogp','dlog WageBar (adj. for schooling)',1,0);
chadfig2('Change in log p','Change in log Earnings (adj. for schooling)',1,0);
makefigwide;
print('-depsc','estimatedeltaYWM_LongDiff');
%title('Long Difference, 1960 - 2010');
print('-dpsc','-append',['estimatedeltaYWM_' CaseName]);
tlestr='Estimating delta*mu from Long Difference, 1960-2010';
[bcoeff,se,tstat,robse,trob]=ols2(dlogwagebarhat,[ones(Noccs,1) dlogp],tlestr,'dlogwagebarhat',strmat('Constant dlogpYWM'));
disp 'Note: the slope is an estimate of delta / (theta(1-eta))';
DeltaOLSLongDiff=bcoeff(2)/mu;
[bcoeff,se,tstat,robse,trob]=iv2(dlogwagebarhat,[ones(Noccs,1) dlogp],[ones(Noccs,1) log(A(:,6)./A(:,1))],tlestr,'dlogwagebarhat',strmat('Constant dlogpYWM'),strmat('Constant dlog(A)'),1);
DeltaIVLongDiff=bcoeff(2)/mu;
DeltaSE=robse(2)/mu;


disp ' '
disp '========================================================';
disp '         SUMMARY OF ESTIMATEDELTAYWM EVIDENCE';
disp '========================================================';
disp ' ';
disp 'Note: the regression slopes above estimate delta / (theta(1-eta))';
disp 'The numbers below are the correct estimates of delta';
disp ' ';
disp '*** OLS: ***';
cshow(' ',[Decades DeltaOLS],'%6.0f %10.4f','Decade delta se tstat');
disp ' ';
dltaOLS=median(DeltaOLS(:,1));
fprintf('The median value of delta across decades is %10.4f\n',dltaOLS);
fprintf('The delta estimtate from the 1960-2010 long diff is %10.4f\n',DeltaOLSLongDiff);
disp ' ';

disp ' ';
disp '*** IV: ***';
cshow(' ',[Decades DeltaIV],'%6.0f %10.4f','Decade delta se tstat');
disp ' ';
dltaIV=median(DeltaIV(:,1));
fprintf('The median value of delta across decades is %10.4f\n',dltaIV);
fprintf('Estimtate from 1960-2010 long diff: %10.4f (%6.4f)\n',[DeltaIVLongDiff DeltaSE]);
disp ' ';

%dlta=dltaIV
dlta=DeltaIVLongDiff

