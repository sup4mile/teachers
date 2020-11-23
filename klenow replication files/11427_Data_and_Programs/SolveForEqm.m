function [YModel,YMktModel,EarningsModel,YwkrModel,LFPModel,ConsumpYoungModel,EarningsYoungModel,GDPYoungModel,WageGapModel,EarningsModel_g,Utility,wModel,HModel,HModelAll,pModel,ExitFlag]=SolveForEqm(TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar);

% SolveForEqm.m    
%
% Given TauW, TauH, Z, etc. this function solves for the general equilibrium.
%
%   Y    = GDP per person = Y in the model
%   Ywkr = GDP per worker = Y / LFP
%   LFP  = Aggregate LFP rate = Fraction of total population that is working
%
%  See 2015-06-02-SolvingGE.pdf notes.
%
%  Method: For each year,
%    1. Guess values for {mgtilde}, Y ==> 5 unknowns
%    2. Solve for {wi} from Hi^supply = Hi^demand
%    3. Compute mghat, Yhat 
%    4. Iterate until converge.

global Noccs Ngroups Ncohorts Nyears CohortConcordance TauW_Orig pData 
global TauW_C phi_C mgtilde_C w_C GeoArith_adjust mgEstimated StopHere % For keeping track of history in solution

StopHere=0;

% Initialize cohort history, needed for solution.
w_C=zeros(Noccs,Ncohorts); 
mgtilde0=11000; % mgtilde := mg^(1/theta*1/(1-eta)) -- better scaled
mgtilde_C=ones(Ngroups,Ncohorts)*mgtilde0; 
TauW_C=zeros(Noccs,Ngroups,Ncohorts); % Cohort order
for g=1:Ngroups;
    TauW_C(:,g,1:6)=flipud(squeeze(TauW(:,g,:))')';
end;


% Guesses
Y0=25000;
%x0=[11000 8000 8000 8000 Y0];
%x0=[1100 800 800 800 Y0]; % For theta(1-eta)=3.44 and eta=1/4
%x0=[1100 800 800 800 Y0]/3; % For theta(1-eta)=3.44 and eta=1/4
%x0=[15000 8000 10000 8000 Y0]; % For theta(1-eta)=1.9 and eta=.10 See mgtilde_C(:,c)=x(1:4); line below for help
%if dlta==1; % Get initial guesses based on mg in EstimateTauZ
    %x0=[[9 1 5 1]*10^7 Y0] % For theta(1-eta)=1.9 and eta=.10 See mgtilde_C(:,c)=x(1:4); line below for help
    %    x0=[[9 1 5 1]*10^6 Y0]; % For theta(1-eta)=1.9 and eta=.10 See mgtilde_C(:,c)=x(1:4); line below for help
    %end;
x0=[(mgEstimated(:,1).^mu)' Y0]; 
%options=optimset('Display','iter'); %,'Algorithm','trust-region-dogleg');
options=optimset('Display','none'); %,'Algorithm','trust-region-dogleg');
wModel=zeros(Noccs,Nyears);
HModel=zeros(Noccs,Nyears);
LFPModel=zeros(Nyears,1);
pModel=zeros(Noccs,Ngroups,3,Nyears); % YMO is 3rd dimension
HModelAll=zeros(Noccs,Ngroups,3,Nyears); % YMO is 3rd dimension
EarningsModel=zeros(Nyears,1);    % Total MARKET Labor Earnings (differs from Y if Revenue~=0)
EarningsModel_g=zeros(Ngroups,Nyears);    % Total MARKET Labor Earnings by Group
%EarningsAllModel=zeros(Nyears,1); % Total Labor Earnings if everyone worked (Pwork=1)
ConsumpYoungModel=zeros(Nyears,1);
EarningsYoungModel=zeros(Nyears,1); % Market (Earnings is always market earnings in model)

% HomeConsumpYModel=zeros(Nyears,1);
% FullConsumpY=zeros(Nyears,1);
% HomeEarningsYModel=zeros(Nyears,1); % Earnings in home sector by the Young
% Home_and_MktOutputY=zeros(Nyears,1); % Home+Mkt Production by the Young
%Home_and_MktConsumpY=zeros(Nyears,1); % Home+Mkt Production by the Young
Utility=zeros(Nyears,1); 
YModel=zeros(Nyears,1); % Includes home occ
YMktModel=zeros(Nyears,1); % Just market occs = sum(w*H) for market occs
GDPYoungModel=zeros(Nyears,1);
WageGapModel=zeros(Ngroups,Nyears); % Just market occs                                    
%WageGapAllModel=zeros(Ngroups,Nyears); % WageGap if everyone worked (Pwork=1)
ExitFlag=zeros(Nyears,1);
xguess=x0;

% Iterate over time to solve the model decade by decade
for t=1:Nyears;
    fprintf('.');
    fune_solve=@(x) e_solveeqm(x,t,TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar);
    [x,fval,flag]=fsolve(fune_solve,xguess,options);
    ExitFlag(t)=flag; if flag~=1; disp 'Exit Flag not equal to one. Stopping...'; keyboard; end;
    [resid,wt,Ht,Yhat,Ymkt,HModelAllt,pModelt,LaborIncome_mkt,LaborIncome_gmkt,WageGapt,ConsumpYoung_t,EarningsYoung_t,GDPYoung_t,Utilityt]=e_solveeqm(x,t,TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar); 
    c=7-t;
    mgtilde_C(:,c)=x(1:4); % Drop the semicolon here to see the mgtilde if x0 guess is wrong
    xguess=x; % Update starting values to previous year's solution
    w_C(:,c)=wt;
    wModel(:,t)=wt;
    HModel(:,t)=Ht;
    HModelAll(:,:,:,t)=HModelAllt;
    pModel(:,:,:,t)=pModelt;
    EarningsModel(t)=LaborIncome_mkt;
    EarningsModel_g(:,t)=LaborIncome_gmkt;
    YModel(t)=Yhat;
    YMktModel(t)=Ymkt;
    WageGapModel(:,t)=WageGapt;
    ConsumpYoungModel(t)=ConsumpYoung_t;
    EarningsYoungModel(t)=EarningsYoung_t;
    GDPYoungModel(t)=GDPYoung_t;
    %HomeEarningsYModel(t)=HomeEarningsYt;
    %Home_and_MktConsumpY(t)=Home_and_MktConsumpYt;
    Utility(t)=Utilityt;

    % Compute LFP and GDP per worker from the solution
    % First, we add up across YMO
    Pwork_gYMO=squeeze(sum(pModelt(2:Noccs,:,:))); % Ngroups x YMO
    Pwork_g=zeros(Ngroups,1);
    %NumPeople_g=zeros(Ngroups,1);
    for ymo=0:2;
        Pwork_g    =Pwork_g+Pwork_gYMO(:,1+ymo).*q(:,c+ymo,t);
        %NumPeople_g=NumPeople_g+q(:,c+ymo,t);
    end;
    LFPModel(t)=sum(Pwork_g); % No need to multiply by NumPeople_g, as we've already done that!
    
    if flag==1;
        x0=x; % Use most recent results for new starting point
    end;
end;
YwkrModel=YMktModel./LFPModel;  % GDP(mkt) per worker = GDP(mkt) per person * Persons/Wkrs


% -------------------------------------------------------
% Sub-functions
% -------------------------------------------------------

% e_solveeqm(x,t,TauH,TauW...)
function [resid,wt,Ht,Yhat,Ymkt,HModelAllt,pmodelt,LaborIncome_mkt,LaborIncome_gmkt,WageGapt,ConsumpYoung_t,EarningsYoung_t,GDPYoung,Utility]=e_solveeqm(x,t,TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar); 

% Given a guess for x=[mgtilde Y] 5x1 and a year t (e.g. 1=1960), 
% solve for w(i) in year t and then compute key moments:
%
%    1. Guess values for {mgtilde}, Y ==> 5 unknowns
%    2. Solve for {wi} from Hi^supply = Hi^demand
%    3. Compute mghat, Yhat 
%    4. Iterate until converge.


global Noccs Ngroups Ncohorts Nyears CohortConcordance TauW_Orig pData ShortNames CaseName
global TauW_C phi_C mgtilde_C w_C pModel % For keeping track of history in solution

mgtilde_t=x(1:4);
Y_t=x(5);
Mkt=2:Noccs;

[wt,Ht,wtildet,HModelAllt,pmodelt]=solveeqm(x,t,TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar); % Returns Noccs x 1 vectors of wages and H(i,t)  and NxG wtilde
mgtildehat=(sum(wtildet.^theta)).^mu;
rho=1-1/sigma;
mu=1/theta*1/(1-eta);
Yhat=sum( (A(:,t).*Ht).^rho ).^(1/rho);
resid=[mgtilde_t-mgtildehat Yhat-Y_t];

% WageGapt (1 x Ngroups) = WageBar(g)/WageBar(WM) net of taxes
%   Average across occupations
%   HModelAllt (Noccs x Ngroups x YMO) but HModelAllt already includes q
%     HAll_i(:,ymo+1)=(q(:,c,t)'.*pig_t.*texp_t.*AvgQuality)';
%   q (Ngroups x Ncohorts x Nyears)
%   pmodelt (Noccs x Ngroups x YMO)  
%   LaborIncomeAll: Total earnings at market prices if *everyone* worked (Pwork=1)
%     From 2016-02-17-EarningsAll.pdf notes, we multiply HModelAllt by Pwork^(mu-1).
%     The mu exponent gets the "per worker" version right, and the -1 adjusts for 
%     the aggregation, converting p(i,g) into ptilde.
%
%   Simple aggregates in the model are all "per person" since our economy has a population=1.

TotalEarnings_ig = mult(squeeze(1-TauW(:,:,t)).*sum(HModelAllt,3),wt);  % Noccs x Ngroups
LaborIncome_gmkt=sum(TotalEarnings_ig(Mkt,:)); % 1 x Ngroups
LaborIncome_mkt=sum(sum(TotalEarnings_ig(Mkt,:)));

% Market GDP = sum(w*H) for mkt occs, adjusting for relative price (2017-11-27)
LaborPayments_ig = mult(sum(HModelAllt,3),wt);  % Noccs x Ngroups
LaborPayments_mkt=sum(sum(LaborPayments_ig(Mkt,:)));
Ymkt = ( LaborPayments_mkt / Y_t^(1/sigma) )^(sigma/(sigma-1));
Pmkt = (Ymkt/Y_t)^(-1/sigma); % Relative price of market goods - for constructing GDPYoung later below

lfp_gc=squeeze(sum(pmodelt(2:Noccs,:,:))); % This is Pwork_gc Ngroups x YMO
c=7-t; co=c+2;
NumWorkers_g = sum(q(:,c:co,t).*lfp_gc,2)';  % 1xG
WageBar_g  = sum(TotalEarnings_ig(Mkt,:))./NumWorkers_g; %1xG  Market!
WageGapt=WageBar_g/WageBar_g(1); % Market

% Consumption: Updated 6/9/16. See Chad-TalentNotes.pdf (page 2c)
% (for all occs including home since a welfare measure)
%     c* = 1/3*(1-eta)*LifetimeIncome = cYoung
%     e*(1+tauh) = eta*LifetimeIncome
EarningsYoung_ig = mult(squeeze(1-TauW(:,:,t)).*HModelAllt(:,:,1),wt);  % Noccs x Ngroups
LIYoung_ig=Tbar(:,:,t)./TExperience(:,:,c,t).*EarningsYoung_ig; % Noccs x Ngroups
NumYoungt=sum(q(:,7-t,t));
cY_ig=1/3*(1-eta)*LIYoung_ig; % Aggregate consumption for the young in i,g since HModelAll aggregates                              
%cY_ig(1,:)=0; % Zero out home values in case they are NaN to make summations work easily.
ConsumpYoung_t=sum(sum(cY_ig))/NumYoungt; % Per young person
if isnan(ConsumpYoung_t); disp 'ConsumpYoung_t isnan'; keyboard; end;
EarningsYoung_t=sum(sum(EarningsYoung_ig(Mkt,:)))/NumYoungt; % Per young person (mkt-based measure)

% GDPYoung (market). See 2017-11-27-MarketGDP.pdf notes
LaborPayYoung_ig = mult(HModelAllt(:,:,1),wt);  % HModelAllt already includes the q(.) and p(.)
LaborPayYoung_mkt =sum(sum(LaborPayYoung_ig(Mkt,:))); 
GDPYoung = LaborPayYoung_mkt/Pmkt / NumYoungt;% Per young Person (market)

% UTILITY: Have to be careful here -- expected utility means an inequality term?
%   In version 4.0, we were ignoring this, so not true utility.
%   Old: We feed in log (c-average) rather than average of log(c) across people. 
%   1/14/19: Now using Chang's formula -- mgtilde! this should be correct!
s_c=1./(1+(1-eta)/beta./phi(:,t));
NumYoung_ig=mult(pmodelt(:,:,1),q(:,c,t)');
Cterm=beta*log(cY_ig./NumYoung_ig);
% Old -- ignores the taste heterogeneity
%Utility_ig = kron(log(1-s_c),ones(1,Ngroups)) + Cterm + log(Z(:,:,c));
%Utility_g2=nansum(NumYoung_ig.*Utility_ig);                 % 1xNgroups

% Utility is log(U) from paper; makes CE-welfare easy to calculate (beta drops out)
Utility_g=log(mgtildehat); % From Chang's 1/14/19 notes -- includes taste heterogeneity as well
NumYoung_g=sum(NumYoung_ig); % 1xG
Utility  =nansum((NumYoung_g.*Utility_g)')/NumYoungt;
Utility=exp(Utility); % To feed the right number in for chaining (like Consumption, not log(c))
