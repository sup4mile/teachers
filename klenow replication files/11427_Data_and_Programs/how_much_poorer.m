function [growthshares,Y_alt,YBaseline,pModel_alt]=how_much_poorer(WhatUnchanged,TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar,GroupToChange,ShowTimeBreakdown);
%function []=how_much_poorer(WhatUnchanged,TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar,GroupToChange,ShowTimeBreakdown);
%
% Details of how_much_poorer.
%
%  WhatUnchanged={'TauWTauH','TauW','TauH','Both+Z','MeansOnly','DispersionOnly'} determines how much is allowed to change
%  GroupToChange={'WW','BM','BW'} to study effect of changing Taus for just that group only.

global Nyears Noccs Ncohorts Ngroups Decades GroupNames ExperienceCohortFactor HAllData CaseName WhatToChain earningsweights_avg ConstrainTauH

if ~exist('GroupToChange'); GroupToChange=0; end;
if ~exist('ShowTimeBreakdown'); ShowTimeBreakdown=0; end;
GetMeans=0;
if isequal(WhatUnchanged,'TauW');         ChangeTauH=0; ChangeTauW=1; ChangeZ=0; end;
if isequal(WhatUnchanged,'TauH');         ChangeTauH=1; ChangeTauW=0; ChangeZ=0; end;
if isequal(WhatUnchanged,'TauWTauH');     ChangeTauH=1; ChangeTauW=1; ChangeZ=0; end;
if isequal(WhatUnchanged,'Both+Z');       ChangeTauH=1; ChangeTauW=1; ChangeZ=1; end;
if isequal(WhatUnchanged,'NoChange');     ChangeTauH=0; ChangeTauW=0; ChangeZ=0; end;
if isequal(WhatUnchanged,'MeansOnly');      ChangeTauH=0; ChangeTauW=0; ChangeZ=0; GetMeans=1; end;
if isequal(WhatUnchanged,'DispersionOnly'); ChangeTauH=0; ChangeTauW=0; ChangeZ=0; GetMeans=1; end;
Mkt=2:Noccs;

disp ' '; disp ' ';
disp '========================================================================';
if GroupToChange==0;
    disp ([' How much poorer in 2010 with 1960 values of ' WhatUnchanged]);
else;
    disp ([' How much poorer in 2010 with 1960 values of ' WhatUnchanged ', for ' GroupNames{GroupToChange}]);
end;    
disp '========================================================================';


% BASELINE -- First we get the baseline values, useful in both sides of how_much_poorer
%YBaseline=SolveForEqm(TauH,TauW,Z,TgHome,TExperience,TigYMO,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar);
load(['SolveEqmBasic_' CaseName]);

pBaseline=pModel; % Baseline model allocation
WageGapBaseline=WageGapBaseline';
WageGapBaseline(:,1)=[];
YBaseline=[GDPBaseline GDPMktBaseline GDPwkrBaseline LFPBaseline EarningsBaseline ConsumpYoungBaseline EarningsYoungBaseline GDPYoungBaseline WageGapBaseline EarningsBaseline_g' UtilityBaseline];
chaintle={'GDP per person','GDP (mkt) per person','GDP (mkt) per worker','Labor Force Participation (LFP)','Earnings (market)','ConsumpYoung (market)','EarningsYoung (market)','GDPYoung (market)','WageGapWW','WageGapBM','WageGapBW','EarningsWM','EarningsWW','EarningsBM','EarningsBW','Utility: CE-Welfare of Young'};

% Means -- for both "MeansOnly" and "DispersionOnly" cases
%  Means of log(1-TauW) and log(1+TauH) keeps it easy
%  N.B.: Only for mkt occs since tauw/tauh in Home are always zero!
if GetMeans; 
    lnTauWmeans=zeros(size(TauW));
    lnTauHmeans=zeros(size(TauH));
    for g=1:Ngroups; % Note well; earningweights_avg is Noccs x 1 -- no time dimension
        lnTWmean_g=nansum(mult(squeeze(log((1-TauW(Mkt,g,:)))),earningsweights_avg(Mkt)));
        lnTHmean_g=nansum(mult(squeeze(log(1+TauH(Mkt,g,:))),earningsweights_avg(Mkt)));
        lnTauWmeans(Mkt,g,:)=kron(lnTWmean_g,ones(Noccs-1,1));
        lnTauHmeans(Mkt,g,:)=kron(lnTHmean_g,ones(Noccs-1,1));
    end;
    lnTauHmeans(:,:,7)=lnTauHmeans(:,:,6); % Fix 7/8 for "prev" analysis below
    lnTauHmeans(:,:,8)=lnTauHmeans(:,:,6);
end;

% For Z's and Tau's we have to loop since they are Noccs x Ngroups x Nyear/cohorts
Z1960=ones(size(Z)); 
TauH1960=ones(size(TauH)); 
TauW1960=ones(size(TauW));
Cohort1960=7-1; % Cohort born in 1960
Cohort2010=7-6; % Cohort born in 2010
for g=1:Ngroups;
    Z1960(:,g,:)    =mult(ones(Noccs,Ncohorts),Z(:,g,Cohort1960));
    TauH1960(:,g,:) =mult(ones(Noccs,Ncohorts),TauH(:,g,Cohort1960));
    TauW1960(:,g,:) =mult(ones(Noccs,Nyears),TauW(:,g,1));
end;    


% Now update the appropriate variables
TauW_alt=TauW;
if ChangeTauW; 
    if GroupToChange==0;
        TauW_alt=TauW1960;
    else % For Group-Specific changes
        TauW_alt(:,GroupToChange,:)=TauW1960(:,GroupToChange,:); 
    end;
end;

TauH_alt=TauH; Z_alt=Z; 
if ChangeTauH; 
    if GroupToChange==0;
        TauH_alt=TauH1960;
    else;
        TauH_alt(:,GroupToChange,:)=TauH1960(:,GroupToChange,:); 
    end;
end;
if ChangeZ; 
    if GroupToChange==0;
        Z_alt=Z1960;
    else;
        Z_alt(:,GroupToChange,:)=Z1960(:,GroupToChange,:); 
    end;
end;  
  
% Dispersion: How much poorer today if 1960 levels of dispersion forever?
% (i.e. let means evolve as they did)
if isequal(WhatUnchanged,'DispersionOnly'); 
    if GroupToChange==0;
        for t=1:Nyears;
            % Careful -- use lnTauWmeans so we need to use log
            ct=7-t;
            TW=log(1-TauW(:,:,1))-lnTauWmeans(:,:,1) + lnTauWmeans(:,:,t); % Update mean only
            TH=log(1+TauH(:,:,Cohort1960))-lnTauHmeans(:,:,Cohort1960) + lnTauHmeans(:,:,ct); 
            TauW_alt(:,:,t)=1-exp(TW);
            TauH_alt(:,:,ct)=exp(TH)-1;
        end;
    else % For Group-Specific changes
        for t=1:Nyears;
            ct=7-t;
            TW=log(1-TauW(:,GroupToChange,1))-lnTauWmeans(:,GroupToChange,1) + lnTauWmeans(:,GroupToChange,t); % Update mean only
            TH=log(1+TauH(:,GroupToChange,Cohort1960))-lnTauHmeans(:,GroupToChange,Cohort1960) + lnTauHmeans(:,GroupToChange,ct); 
            TauW_alt(:,GroupToChange,t)=1-exp(TW);
            TauH_alt(:,GroupToChange,ct)=exp(TH)-1;
        end;
    end;      
    TauH_alt(TauH_alt<ConstrainTauH)=ConstrainTauH; % Enforce the constraint on counterfactuals as well
    TauW_alt(TauW_alt>0.99)=0.99; % Cannot have >=1
end;


% Means: How much poorer today if 1960 mean tau's forever?
% (i.e. let dispersion evolve as they did)
if isequal(WhatUnchanged,'MeansOnly'); 
    if GroupToChange==0;
        for t=1:Nyears;
            % Careful -- use lnTauWmeans so we need to use log
            ct=7-t; % No need to worry about cohorts 7 and 8 since only care about 2010
            TW=log(1-TauW(:,:,t))-lnTauWmeans(:,:,t) + lnTauWmeans(:,:,1); % means unchanged
            TH=log(1+TauH(:,:,ct))-lnTauHmeans(:,:,ct) + lnTauHmeans(:,:,Cohort1960); 
            TauW_alt(:,:,t)=1-exp(TW);
            TauH_alt(:,:,ct)=exp(TH)-1;
        end;
    else % For Group-Specific changes
        for t=1:Nyears;
            ct=7-t; % No need to worry about cohorts 7 and 8 since only care about 2010
            TW=log(1-TauW(:,GroupToChange,t))-lnTauWmeans(:,GroupToChange,t) + lnTauWmeans(:,GroupToChange,1); % 1960 means
            TH=log(1+TauH(:,GroupToChange,ct))-lnTauHmeans(:,GroupToChange,ct) + lnTauHmeans(:,GroupToChange,Cohort1960); 
            TauW_alt(:,GroupToChange,t)=1-exp(TW);
            TauH_alt(:,GroupToChange,ct)=exp(TH)-1;
        end;
    end;      
    % Note that this changes the variance slightly!
    TauH_alt(TauH_alt<ConstrainTauH)=ConstrainTauH; % Enforce the constraint on counterfactuals as well
    TauW_alt(TauW_alt>0.99)=0.99; % Cannot have >=1
end;

%checktaudispersion_alt


t=6;
[y_output,y_mkt,y_earnings,y_wkr,lfp,consumpmkt,earningsyoung,gdpyoung,wagegap,earnings_g,utility,wModel,HModel,HModelAll,pModel]...
    =SolveForEqm(TauH_alt,TauW_alt,Z_alt,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar);
Yinit=[y_output(t) y_mkt(t) y_wkr(t) lfp(t) y_earnings(t) consumpmkt(t) earningsyoung(t) gdpyoung(t) wagegap(2:4,t)' earnings_g(:,t)' utility(t)];

gIactual=log(YBaseline(t,:)./YBaseline(1,:))/(Decades(end)-Decades(1));
gImodel =log(Yinit./YBaseline(1,:))/(Decades(end)-Decades(1));
gIshare = (1-gImodel./gIactual)*100;
rowtle=['-- 1960 --';         '2010(alt) ';         '-- 2010 --'];
fmt='%8.0f %8.0f %8.0f %8.3f %8.0f %8.0f %8.0f %8.0f %8.3f %8.3f %8.3f %8.0f %8.0f %8.0f %8.0f %8.0f';

cshow(rowtle,[YBaseline(1,:); Yinit; YBaseline(t,:)],fmt,'Y Ymkt Ywkr LFP Earn Cons EarnY gdpY GapWW GapBM GapBW EarnWM EarnWW EarnBM EarnBW Util');
cshow('Growth(alt)',gImodel,'%8.4f',[],'nonee',1);
cshow('Growth(tru)',gIactual,'%8.4f',[],'nonee',1);
cshow('Difference ',gIactual-gImodel,'%8.4f',[],'nonee',1);
cshow('// Share //',gIshare,'%8.1f',[],'nonee',1);
disp ' '; disp ' ';

growthshares=gIshare; % To return from the function
pModel_alt=pModel;
Y_alt = [y_output y_mkt y_wkr lfp y_earnings consumpmkt earningsyoung gdpyoung wagegap(2:4,:)' earnings_g'];
%           1       2     3    4     5          6            7           8       9:11            12:15
%keyboard


