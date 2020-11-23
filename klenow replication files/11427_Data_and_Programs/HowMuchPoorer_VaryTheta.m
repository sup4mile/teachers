% HowMuchPoorer_VaryTheta.m   
%
%  HowMuchPoorer basic approach:
%  Alternative approach to evaluating contribution of tau's etc.
% "Relative to 2010 in model solution, how much poorer would we be
%  if ______ had occurred?"
%
%  Do this with the *benchmark* tau's for different values of THETA
%  -- i.e. do not re-estimate the tau's as we change theta.
%
%  How does the growth contribution of tauw/tauh look in that case?

clear; global CaseName;
diarychad('HowMuchPoorer_VaryTheta');
CaseName='Benchmark' % Use the benchmark A's, phi's, tau's

global Noccs Ngroups Ncohorts Nyears CohortConcordance TauW_Orig pData HAllData Decades ExperienceCohortFactor
global TauW_C phi_C mgtilde_C w_C WhatToChain % For keeping track of history in solution

load(['TalentData_' CaseName]); % From EstimateTauZ and earlier programs

ThetaValues=[1.5 1.75 1.9 2 2.25 2.5 3]'% 4]'
ChainSingleCase=1
GrowthShare=zeros(length(ThetaValues),1);

for i=1:length(ThetaValues);
    theta=ThetaValues(i);
    ShowParameters;
    [GrowthShare_TWTH]=how_much_poorer('TauWTauH',TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar);
    GrowthShare(i)=GrowthShare_TWTH(1);
end;

cshow(' ',[ThetaValues GrowthShare],'%12.3f','Theta GrowthShare');

definecolors;
figure(1); figsetup;
plot(ThetaValues,GrowthShare,'-','Color',myblue);
hold on;
plot(ThetaValues,GrowthShare,'o','Color',mygreen);
chadfig2('Theta','Share of growth due to tau_w and tau_h',1,0)
print HowMuchPoorer_VaryTheta.eps
    
save(['HowMuchPoorer_VaryTheta']);
diary off;


