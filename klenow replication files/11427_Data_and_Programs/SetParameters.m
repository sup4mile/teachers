% SetParameters.m
%
%  Sets the Benchmark basic parameters for the Talent project

% Make sure all relevant files are in path
curpath=path;
if isempty(strfind(curpath,'ChadMatlab'));% | isempty(strfind(curpath,'Work/procs'));
    if exist('ChadMatlab')==7; 
        curdir=pwd;
        path([curdir '/ChadMatlab'],path); 
    end; 
end;

% Parameters:

    dlta=0;    % Share of "taste" types
    eta=0.103;
    theta = 2; % New benchmark 12/15/17 %  1.517; % Via Pete 11/27/17
    %theta=KeyThetaParam/(1-eta);
    
    % On beta: Weight in utility on log(c) versus log(1-s) - 3 period model
    % In the new 4/4/16 lifecycle version, beta is replaced by 3*beta, so I'm just
    % making that adjustment here rather than replacing beta with 3*beta in all programs
    % that follow ==> it stays at 0.693.
    beta=3*0.693/3;  
    sigma=3;
 
    HOME=1;
    WhatToChain='Output'; % 'Output' and 'Earnings' are the alternatives
    FiftyFiftyTauHat=0; % Turn this on for robustness to split tauhat 50/50 in every year
    IgnoreBrawnyOccupations=0; % Turn this on for robustness to zero out tauh/w in brawny occupations
    NoFrictions2010=0;  % Turn this on for robustness to choose T(i,g,2010) s.t. set tauw(2010)=tauh(2010)=0.    
    WageGapAdjustmentFactor=1; % Fraction of Wage Gap to preserve. Zero ==> use Earnings(WM)
                               
    % Nov 2017: New benchmark case has fixed split of 50/50 for 1960, no iterating
    %Alpha0FixedSplit=NaN;  % Alternative is e.g. AlphaFixedSplit=0.5 to use that AlphaSplitTauW1960 value and not iterate
    Alpha0FixedSplit=0.5; % 1/2 split initially
    AlphaSplitTauW1960=0.5; % Default starting value
    ConstantAlpha = NaN; % For implementing fixed splits other than 50/50
    NoGroupExpAdj=0;      % Default is to adjust wage growth for Y->M by GroupAdjExperience
    WhichWageGrowth=0;    % Default is to use both Y->M and M->O wage growth when estimating TauW (1=Y->M, 2=M->O)
    HalfExperience=0;     % Default is WM and WW get same return to experience. Turning this on gives 1/2 return to WW/BM/BW
    
    OccupationtoIdentifyWageHome=23; % WageBarHome(WM)=WageBar(Secretaries,WM)
    
    ChainSingleCase=0; % Turn on if we only wish to chain the TauWTauH case
    HighQualityFigures=0; % Turn on to place labels in certain figures ==> programs wait for user input
    SameExperience=1;  % Defaults is same returns to experience in all occs
    ConstantExperience=0; % If 1, make the returns to experience constant over time
    ConstrainTauH=-0.8;  % Lower bound on how negative TauH can get (meaningless if below -1). -999=unconstrained
    NumHomeDraws=[]; % Legacy parameter, not used.
    PurgeWageGrowthSelection=0; % Default is *not* to do anything here; seems perverse
    EstimateDelta=0; % Turn on to estimate delta using YWM wagehat vs p