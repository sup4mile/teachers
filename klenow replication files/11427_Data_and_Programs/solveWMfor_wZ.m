
% solveWMfor_wZ.m: Estimate w(i) and ZWM(i) 
%   Given phi/s from EducationYWM.
%   Requires TWO normalizations (because we have N-1 wages since no wage observed at home
%   and N+1 unknowns (mwm and Z in every occupation)
%     1. Zhome=1
%     2. wagebarHome=wagebar(OccupationtoIdentifyWageHome)


function [W,PHI,S,MWM,ZWM,WageBarHomeY]=solveWMfor_wZ(eta,beta,theta,gam,sigma,OccupationtoIdentifyWageHome,Tbar,ShowData);

global CaseName;
load(['CohortData_' CaseName '.mat']);

HOME=1;
Mkt=2:Noccs; % Market Occupation numbers

% Requires Ty/Tbar in order to use the wagebar equation.
Ty=ones(Noccs,1); % Implicitly assumed below, not passed to the functions.
Tbar=squeeze(Tbar(:,1,:)); % Pull WM ==> Noccs x Nyears

% Store results in Noccs x Decades matrices
W=zeros(Noccs,Nyears)*NaN;
MWM=zeros(Nyears,1)*NaN;   % Nyears x 1
WageBarHomeY=zeros(Ngroups,Nyears)*NaN;
mu=1/theta*1/(1-eta);
etabar=eta^(eta/(1-eta));

% Use Education to get phi and s
S=EducationYWM/25; % 25 years is denominator
PHI=(1-eta)/beta*S./(1-S); % Noccs x Nyears



for t=1:Nyears;

    disp ' '; disp ' ';
    disp '-------------------------------------------------------------';
    fprintf('Estimating w and Z from YWM: %6.0f\n',CohortConcordance(t,1));
    disp '-------------------------------------------------------------';

    YoungCohort=CohortConcordance(t,2);
    wagebar=squeeze(WageBar(:,:,YoungCohort,t)); 
    pt=squeeze(p(:,:,YoungCohort,t));
    phi=PHI(:,t);
    s=S(:,t);

    % First get wagebar(HOME,WM) and then Zhome=1 to get MWM
    wagebar(HOME,WM)=wagebar(OccupationtoIdentifyWageHome,WM);
    wagebarHomeYWM=wagebar(HOME,WM);
    mwm = pt(HOME,WM)^(-dlta)*(wagebarHomeYWM*(1-s(HOME))^(1/beta) / (gam*etabar*Ty(HOME)/Tbar(HOME,t)) ).^(1/mu);

    % Now, Z comes from the wagebar equation for each occupation.
    z=( gam*etabar*(pt(:,WM).^dlta*mwm).^mu.*Ty./Tbar(:,t)./wagebar(:,WM) ).^beta ./(1-s);

    % And w from the wtilde expression 
    wtilde = (mwm*pt(:,WM)).^(1/theta); % Noccs x 1
    w = wtilde ./ (Tbar(:,t).*s.^phi.*((1-s).*z).^((1-eta)/beta)); % Noccs x 1
        
    
    fmt='%8.3f %8.3f %10.1f %8.3f %8.3f %8.3f %8.1f %8.1f %8.3f';
    tle='s phi w zWM';
    cshow(OccupationNames,[s phi w z],fmt,tle);
    fprintf('\n  MWM = %12.3f\n',mwm);
    
    % Wagebar in Home for other groups (since Zhome=1 for all groups)
    wagebarHome_g =wagebarHomeYWM*(pt(HOME,:)/pt(HOME,WM)).^(-(1-dlta)*mu);
    
    W(:,t)=w; 
    ZWM(:,t)=z;
    MWM(t)=mwm; 
    WageBarHomeY(:,t)=wagebarHome_g';
    
end; % Looping over decades

disp ' '; disp ' ';
disp '=============================================================================';
disp ' Summary of Results for w, phi, s, Z from Young WM... (solveWMfor_w_phi.m)'
disp '=============================================================================';
disp ' ';
disp 'Values for w(i,t) --- wage per unit of human capital';
tle='1960 1970 1980 1990 2000 2010';
cshow(OccupationNames,W,'%8.0f',tle);

disp ' '; disp ' ';
disp 'Values for phi(i,t)';
tle2='1960 1970 1980 1990 2000 2010';
cshow(OccupationNames,PHI,'%8.3f',tle2);

disp ' '; disp ' ';
disp 'Values for s(i,t)';
tle2='1960 1970 1980 1990 2000 2010';
cshow(OccupationNames,S,'%8.3f',tle2);


disp ' '; disp ' ';
disp 'Values for ZWM(i,t)';
tle2='1960 1970 1980 1990 2000 2010';
cshow(OccupationNames,ZWM,'%8.3f',tle2);


disp ' '; disp ' ';
disp 'WageBarHomeY(g,t):';
cshow(GroupNames,WageBarHomeY,'%8.0f',tle2);
disp ' '; disp ' ';

