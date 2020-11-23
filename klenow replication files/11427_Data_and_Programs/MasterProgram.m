% MasterProgram.m     "The Allocation of Talent and U.S. Economic Growth"
%
%  Main program to generate matlab results for the Talent paper.
%  Calls the other programs in the correct order and conducts robustness checks
%
%  Note: These cases can be run in parallel, e.g. by starting separate matlab sessions
%        and copying the code into a separate program (or just pasting it in interactively).
%
%  Each case takes about 15 minutes to run. Whole thing is about 4 hours.
%  See the *.log files for results and *.ps and *.eps for figures. 

% % Sandbox
% clear all; global CaseName;
% CaseName='Sandbox';
% SetParameters;
% HighQualityFigures=0;
% % Using Earnings instead of WageBar for GetTExperience
% ReadCohortData
% EstimateTauZ
% CleanandShowTauAZ
% SolveEqmBasic
% HowMuchPoorer

% abc

% % Need to run Additional Figures
% % Benchmark
tic
clear all; global CaseName;
CaseName='Benchmark';
SetParameters;
HighQualityFigures=0;
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
%AdditionalFigures
chadtimer


% EstimateDelta
clear all; global CaseName;
CaseName='EstimateDelta';
SetParameters;
HighQualityFigures=0;
EstimateDelta=1
dlta=0
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer


% Delta=0 
clear all; global CaseName;
CaseName='Delta0';
SetParameters;
HighQualityFigures=0;
dlta=0
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
% AdditionalFigures
%chadtimer







% Delta1 
clear all; global CaseName;
CaseName='Delta1';
SetParameters;
HighQualityFigures=0;
dlta=1
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
% AdditionalFigures
%chadtimer



% DeltaHalf 
clear all; global CaseName;
CaseName='DeltaHalf';
SetParameters;
HighQualityFigures=0;
dlta=1/2
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
% AdditionalFigures
%chadtimer



% HalfExperience: WW/BM/BW get 1/2 the returns to experience of WM Tbar differs by groups
clear all; global CaseName;
CaseName='HalfExperience';
SetParameters;
HalfExperience=1;
HighQualityFigures=0;
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer


% WageGapAdjustmentFactor=1/2
clear all; global CaseName;
CaseName='WageGapHalf';
SetParameters;
WageGapAdjustmentFactor=1/2;
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer

% WageGapAdjustmentFactor=0
clear all; global CaseName;
CaseName='WageGapZero';
SetParameters;
WageGapAdjustmentFactor=0;
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer


% Robustness Check: zero out tauh/w in brawny occupations
clear all; global CaseName;
CaseName='NoBrawny';
SetParameters;
IgnoreBrawnyOccupations=1;
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
chadtimer


% Robustness Check: choosing T(i,g) s.t. set tauw(2010)=tauh(2010)=0
clear all; global CaseName;
CaseName='NoFrictions2010';
SetParameters;
NoFrictions2010 = 1; 
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
chadtimer



% NoConstrainTauH
clear all; global CaseName;
CaseName='NoConstrainTauH';
SetParameters;
ConstrainTauH=-999
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
chadtimer


%  eta=.20
clear all; global CaseName;
CaseName='Eta20';
SetParameters;
eta=.20
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
chadtimer


%  eta=.05
clear all; global CaseName;
CaseName='Eta05';
SetParameters;
eta=.05
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
chadtimer



% Robustness: Theta = 4
clear all; global CaseName;
CaseName='Theta4';
SetParameters;
theta=4;
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
chadtimer

% Robustness: Theta = 1.5
clear all; global CaseName;
CaseName='ThetaLow';
SetParameters;
theta=1.5;
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
chadtimer


% Sigma = 10
clear all; global CaseName;
CaseName='Sigma10';
SetParameters;
HighQualityFigures=0;
sigma=10;
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
chadtimer

% Sigma = 1.05
clear all; global CaseName;
CaseName='SigmaCobb';
SetParameters;
HighQualityFigures=0;
sigma=1.05
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
chadtimer



% AllTauH
clear all; global CaseName;
CaseName='AllTauH';
SetParameters;
ConstantAlpha=0
ConstrainTauH=-2
HighQualityFigures=0;
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
% % AdditionalFigures


% AllTauHconstrained
clear all; global CaseName;
CaseName='AllTauHconstrained';
SetParameters;
ConstantAlpha=0
ConstrainTauH=-0.80
HighQualityFigures=0;
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
% % AdditionalFigures




% AllTauW
clear all; global CaseName;
CaseName='AllTauW';
SetParameters;
ConstantAlpha=1
HighQualityFigures=0;
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
% % AdditionalFigures



% Robustness Check: Splitting TauHat 50/50 in every year into tauw vs tauh
% using ConstantAlpha=1/2

clear all; global CaseName;
CaseName='ConstantAlphaHalf';
SetParameters;
ConstantAlpha=1/2;
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
chadtimer




% Robustness Check: Splitting TauHat 50/50 in every year into tauw vs tauh

clear all; global CaseName;
CaseName='FiftyFifty';
SetParameters;
FiftyFiftyTauHat=1;
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
chadtimer



% WhichWageGrowth=2 (M->O only)
clear all; global CaseName;
CaseName='WageGrowthMO';
SetParameters;
WhichWageGrowth=2; % (M->O only)
Alpha0FixedSplit=0.5; % 1/2 split initially
HighQualityFigures=0;
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
% % AdditionalFigures


% WhichWageGrowth=1 (Y->M only)
clear all; global CaseName;
CaseName='WageGrowthYM';
SetParameters;
WhichWageGrowth=1; % (Y->M only)
Alpha0FixedSplit=0.5; % 1/2 split initially
HighQualityFigures=0;
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
% % AdditionalFigures


% NoGroupExpAdj
clear all; global CaseName;
CaseName='NoGroupExpAdj';
SetParameters;
NoGroupExpAdj=1;
HighQualityFigures=0;
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer


% Alpha0FixedSplit=1/3
clear all; global CaseName;
CaseName='Alpha0is33';
SetParameters;
Alpha0FixedSplit=1/3; % 1/4 split initially
HighQualityFigures=0;
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
% % AdditionalFigures

% Alpha0FixedSplit=2/3
clear all; global CaseName;
CaseName='Alpha0is67';
SetParameters;
Alpha0FixedSplit=2/3; % 
HighQualityFigures=0;
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer


% ConstrainTauH = -0.50 -- try with a *tighter* constraint?
clear all; global CaseName;
CaseName='ConstrainTauHHalf';
SetParameters;
ConstrainTauH=-0.5
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer
chadtimer


% HomeSales
clear all; global CaseName;
CaseName='HomeSales';
SetParameters;
HighQualityFigures=0;
OccupationtoIdentifyWageHome=22; % Sales (22) insteade of default Secretaries (23)
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer



% Delta25 
clear all; global CaseName;
CaseName='Delta25';
SetParameters;
HighQualityFigures=0;
dlta=0.25
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer



% Delta75 
clear all; global CaseName;
CaseName='Delta75';
SetParameters;
HighQualityFigures=0;
dlta=0.75
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer

% Delta90
clear all; global CaseName;
CaseName='Delta90';
SetParameters;
HighQualityFigures=0;
dlta=0.90
ReadCohortData
EstimateTauZ
CleanandShowTauAZ
SolveEqmBasic
HowMuchPoorer



% ChainingVaryTheta
%  - Benchmark parameters/A's/phi's/tau's but vary theta...
HowMuchPoorer_VaryTheta


% % TauWWisZero
% clear all; global CaseName;
% CaseName='TauWWisZero';
% SetParameters;
% ConstantAlpha=1/2; % Try imposing alpha=1/2 in all decades
% ConstrainTauH=-999; % Do not constrain TauH
% HighQualityFigures=0;
% ReadCohortDataWW
% EstimateTauZ
% CleanandShowTauAZ
% SolveEqmBasic
% HowMuchPoorer


% % Trim 1%
% clear all; global CaseName;
% CaseName='Trim1';
% SetParameters;
% HighQualityFigures=0;
% dlta=0
% %phiHome=zeros(1,6); % So that phiHome=phi(23) Secretary
% OccupationtoIdentifyPhi=22; % Sales (22)
% ReadCohortData
% EstimateTauZ
% CleanandShowTauAZ
% SolveEqmBasic
% HowMuchPoorer




