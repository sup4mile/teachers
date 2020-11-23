% HowMuchPoorer.m
%
%  Alternative approach to evaluating contribution of tau's etc.
% "Relative to 2010 in model solution, how much poorer would we be
%  if ______ had occurred?"
%   Ex: if tau's stayed at 1960 level?
%       if tauw had not changed since 1960
%       if tauh had not changed since 1960
%       if only mean tau had changed since 1960, not dispersion
%       if dispersion changed, but not the means
%
%  That is, allow everything else to change as it did, other than ______
%  "// Share //" row reports the contribution in terms of growth rates.
%
clear; global CaseName;
diarychad('HowMuchPoorer',CaseName);

global Noccs Ngroups Ncohorts Nyears GroupNames CohortConcordance TauW_Orig pData HAllData Decades ExperienceCohortFactor
global TauW_C phi_C mgtilde_C w_C WhatToChain earningsweights_avg ConstrainTauH % For keeping track of history in solution

load(['TalentData_' CaseName]); % From EstimateTauZ and earlier programs
ShowParameters;
if exist('NumHomeDraws')~=1; NumHomeDraws=[]; end; %10^5; end;
if isequal(CaseName,'Benchmark');
    ShowTimeBreakdown=1;
else;
    ShowTimeBreakdown=0;
end;

tle={'GDP per person','GDP (mkt) per person','GDP (mkt) per worker','Labor Force Participation (LFP)','Earnings (market)','ConsumpYoung (total)','EarningsYoung (market)','GDPYoung (market)','WageGapWW','WageGapBM','WageGapBW','EarningsWM','EarningsWW','EarningsBM','EarningsBW','Utility: CE-Welfare of Young'};
shorttle=strmat('Y Ymkt Ywkr LFP Earn Cons EarnY gdpY GapWW GapBM GapBW EarnWM EarnWW EarnBM EarnBW Util');
disp 'Concordance of short column titles and concepts:';
for i=1:length(tle);
   disp(['  ' shorttle(i,:) '  ' tle{i}]) 
end;



%how_much_poorer('NoChange',TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar,0,ShowTimeBreakdown);
[growthshareTWTH,Y_TWTH,YBaseline,pModel_TWTH]=how_much_poorer('TauWTauH',TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar,0,ShowTimeBreakdown);

if ~ChainSingleCase;
    how_much_poorer('TauH',TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar,0,ShowTimeBreakdown);
    how_much_poorer('TauW',TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar,0,ShowTimeBreakdown);
    how_much_poorer('Both+Z',TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar,0,ShowTimeBreakdown);
    how_much_poorer('MeansOnly',TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar,0,ShowTimeBreakdown);
    how_much_poorer('DispersionOnly',TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar,0,ShowTimeBreakdown);
end;

% Separate tauw/tauh contributions by Group
if isequal(CaseName,'Benchmark');
    disp ' ';
    disp '***************  BREAKDOWN BY GROUP AND TIME PERIOD *****************'; disp ' ';
    ShowTimeBreakdown=0; % This is not implemented currently
    how_much_poorer('TauWTauH',TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar,WW,ShowTimeBreakdown);
    how_much_poorer('TauH',    TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar,WW,ShowTimeBreakdown);
    how_much_poorer('TauW',    TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar,WW,ShowTimeBreakdown);
    how_much_poorer('TauWTauH',TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar,BM,ShowTimeBreakdown);
    how_much_poorer('TauH',    TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar,BM,ShowTimeBreakdown);
    how_much_poorer('TauW',    TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar,BM,ShowTimeBreakdown);
    how_much_poorer('TauWTauH',TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar,BW,ShowTimeBreakdown);
    how_much_poorer('TauH',    TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar,BW,ShowTimeBreakdown);
    how_much_poorer('TauW',    TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar,BW,ShowTimeBreakdown);
end;
    
save(['HowMuchPoorer_' CaseName]);

diary off;

