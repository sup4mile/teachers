% Variance_IncomeWages.m
%
%  Create data plots of within-occupation variance of ln(Income) or ln(Wages)
%  plotted against relative propensities. In 1980 and Change(1960,2010)

help('Variance_IncomeWages');

disp ' ';
disp ' ';
disp '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
disp '  White Men: Earnings Growth versus Propensity Growth (long difference)'
disp '   11/13/18 -- for Changs shorthand estimation of "alpha" (fraction ability/taste)'
disp '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

% 1960-2010
t60=find(Decades==1960);
cYMO60=CohortConcordance(t60,ymo);
t10=find(Decades==2010);
cYMO10=CohortConcordance(t10,ymo);

pgrowth=log(p(:,WM,cYMO10,t10)./p(:,WM,cYMO60,t60));
wgrowth=log(WageBar(:,WM,cYMO10,t10)./WageBar(:,WM,cYMO60,t60));


figure(2); figsetup;
plotnamesym(pgrowth,wgrowth,ShortNames,9,[],.06,.00);
chadfig2('Change in log pWM, 1960-2010','Change in log wageWM, 1960-2010',1,0);
makefigwide;
print WMWageGrowth1960_2010.eps
cshow(ShortNames,[pgrowth wgrowth],'%8.4f','pGrowth WageGrowth');
ols(wgrowth,[ones(Noccs,1) pgrowth],'ols: WM Wage Growth on pgrowth, 1960-2010','DlnEarnings',strmat('Constant Dlog(p)'));


% 1980-2010
t80=find(Decades==1980);
cYMO80=CohortConcordance(t80,ymo);
t10=find(Decades==2010);
cYMO10=CohortConcordance(t10,ymo);

pgrowth=log(p(:,WM,cYMO10,t10)./p(:,WM,cYMO80,t80));
wgrowth=log(WageBar(:,WM,cYMO10,t10)./WageBar(:,WM,cYMO80,t80));


figure(2); figsetup;
plotnamesym(pgrowth,wgrowth,ShortNames,9,[],.06,.00);
chadfig2('Change in log pWM, 1980-2010','Change in log wageWM, 1980-2010',1,0);
makefigwide;
print WMWageGrowth1980_2010.eps
cshow(ShortNames,[pgrowth wgrowth],'%8.4f','pGrowth WageGrowth');
ols(wgrowth,[ones(Noccs,1) pgrowth],'ols: WM Wage Growth on pgrowth, 1980-2010','DlnEarnings',strmat('Constant Dlog(p)'));


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Var_lnIncome and Var_lnWage versus Relative Propensity Graphs
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%
% INCOME 
% %%%%%%%%%%%%%%%%%

disp ' ';
disp '**********************************************';
disp '          Variance Results: Income';
disp '**********************************************';
disp ' ';

% Var_lnIncome for white women in 1980 vs relative p's
ymo=1+1; % Young
t80=find(Decades==1980);
cYMO80=CohortConcordance(t80,ymo);
pWW80=p(:,WW,cYMO80,t80);
Var_lnInc_WW80=Var_lnIncome(:,WW,cYMO80,t80);
Var_lnInc_WM80=Var_lnIncome(:,WM,cYMO80,t80);
dVar_lnInc80=Var_lnInc_WW80-Var_lnInc_WM80;


% Levels for WW
figure(1); figsetup;
plotnamesym(log(pWW80),Var_lnInc_WW80,ShortNames,9,[],.06,.00);
chadfig2('Absolute propensity, p(ww)','Variance of ln(Income)',1,0);
makefigwide;
xt=[1/64 1/16 1/4 1 4 16 64];
set(gca,'XTick',log(xt));
set(gca,'XTickLabel',strmat('1/64#1/16# 1/4#  1 #  4 # 16 # 64 ','#'));
print VarianceInc_AbsP1980.eps
disp 'Variance of ln(Income) in 1980: WW';
cshow(ShortNames,[Var_lnInc_WW80 pWW80],'%10.4f','Var(WW) p(WW)');
ols(Var_lnInc_WW80,[ones(Noccs,1) log(pWW80)],'ols: Levels (Variance of Income)','Var_lnInc_WW80',strmat('Constant logAbsProp'));


% Diff WW-WM
figure(1); figsetup;
plotnamesym(log(relpWW80),dVar_lnInc80,ShortNames,9,[],.06,.00);
chadfig2('Relative propensity, p(ww)/p(wm)','Difference in Variance of ln(Income) (WW-WM)',1,0);
makefigwide;
xt=[1/64 1/16 1/4 1 4 16 64];
set(gca,'XTick',log(xt));
set(gca,'XTickLabel',strmat('1/64#1/16# 1/4#  1 #  4 # 16 # 64 ','#'));
print VarianceInc_RelP1980.eps
disp 'Variance of ln(Income) in 1980: WW vs WM';
cshow(ShortNames,[Var_lnInc_WW80 Var_lnInc_WM80 dVar_lnInc80 relpWW80],'%10.4f','Var(WW) Var(WM) dVar RelP');
ols(dVar_lnInc80,[ones(Noccs,1) log(relpWW80)],'ols: Levels (Variance of Income)','dVar_lnInc80',strmat('Constant logRelProp'));


% Now the change for Figure 2
t60=find(Decades==1960);
cYMO60=CohortConcordance(t60,ymo);
Var_lnInc_WW60=Var_lnIncome(:,WW,cYMO60,t60);
Var_lnInc_WM60=Var_lnIncome(:,WM,cYMO60,t60);
dVar_lnInc60=Var_lnInc_WW60-Var_lnInc_WM60;

t10=find(Decades==2010);
cYMO10=CohortConcordance(t10,ymo);
Var_lnInc_WW10=Var_lnIncome(:,WW,cYMO10,t10);
Var_lnInc_WM10=Var_lnIncome(:,WM,cYMO10,t10);
dVar_lnInc10=Var_lnInc_WW10-Var_lnInc_WM10;

change_dvar=dVar_lnInc10-dVar_lnInc60;
changelogp=log(relpWW10)-log(relpWW60);

figure(2); figsetup;
plotnamesym(changelogp,change_dvar,ShortNames,9,[],.06,.00);
chadfig2('Change in log p(ww)/p(wm), 1960-2010','Change in Differential Variance(Income), 1960-2010',1,0);
makefigwide;
%xt=[.01 .04 .16 .5 1 4 16 64];
%xt=[1/64 1/16 1/4 1 4 16 64];
%set(gca,'XTick',log(xt));
%set(gca,'XTickLabel',strmat('.01 .04 .16 .5 1 4 16 64'));
%set(gca,'XTickLabel',strmat('1/64#1/16# 1/4#  1 #  4 # 16 # 64 ','#'));
print VarianceInc_RelPChange.eps
cshow(ShortNames,[changelogp change_dvar],'%8.4f','DlnRelP Change_dVar');
ols(change_dvar,[ones(Noccs,1) changelogp],'ols: Changes in Diff Variance of Income','Chang_dVar(Income)',strmat('Constant DlogRelProp'));


% %%%%%%%%%%%%%%%%%
% WAGES 
% %%%%%%%%%%%%%%%%%

disp ' ';
disp '**********************************************';
disp '          Variance results: WAGES';
disp '**********************************************';
disp ' ';

% Var_lnWage for white women in 1980 vs relative p's
ymo=1+1; % Young
t80=find(Decades==1980);
cYMO80=CohortConcordance(t80,ymo);
Var_lnWage_WW80=Var_lnWage(:,WW,cYMO80,t80);
Var_lnWage_WM80=Var_lnWage(:,WM,cYMO80,t80);
dVar_lnWage80=Var_lnWage_WW80-Var_lnWage_WM80;


% Levels for WW
figure(1); figsetup;
plotnamesym(log(pWW80),Var_lnWage_WW80,ShortNames,9,[],.06,.00);
chadfig2('Absolute propensity, p(ww)','Variance of ln(Wage)',1,0);
makefigwide;
xt=[1/64 1/16 1/4 1 4 16 64];
set(gca,'XTick',log(xt));
set(gca,'XTickLabel',strmat('1/64#1/16# 1/4#  1 #  4 # 16 # 64 ','#'));
print VarianceWage_AbsP1980.eps
disp 'Variance of ln(Wage) in 1980: WW';
cshow(ShortNames,[Var_lnWage_WW80 pWW80],'%10.4f','Var(WW) p(WW)');
ols(Var_lnWage_WW80,[ones(Noccs,1) log(pWW80)],'ols: Levels (Variance of Wage)','Var_lnWage_WW80',strmat('Constant logAbsProp'));



figure(1); figsetup;
plotnamesym(log(relpWW80),dVar_lnWage80,ShortNames,9,[],.06,.00);
chadfig2('Relative propensity, p(ww)/p(wm)','Difference in Variance of ln(Wage) (WW-WM)',1,0);
makefigwide;
xt=[1/64 1/16 1/4 1 4 16 64];
set(gca,'XTick',log(xt));
set(gca,'XTickLabel',strmat('1/64#1/16# 1/4#  1 #  4 # 16 # 64 ','#'));
print VarianceWage_RelP1980.eps
disp 'Variance of ln(Wage) in 1980: WW vs WM';
cshow(ShortNames,[Var_lnWage_WW80 Var_lnWage_WM80 dVar_lnWage80 relpWW80],'%10.4f','Var(WW) Var(WM) dVar RelP');
ols(dVar_lnWage80,[ones(Noccs,1) log(relpWW80)],'ols: Levels (Variance of Wage)','dVar_lnWage80',strmat('Constant logRelProp'));


% Now the change for Figure 2
t60=find(Decades==1960);
cYMO60=CohortConcordance(t60,ymo);
Var_lnWage_WW60=Var_lnWage(:,WW,cYMO60,t60);
Var_lnWage_WM60=Var_lnWage(:,WM,cYMO60,t60);
dVar_lnWage60=Var_lnWage_WW60-Var_lnWage_WM60;

t10=find(Decades==2010);
cYMO10=CohortConcordance(t10,ymo);
Var_lnWage_WW10=Var_lnWage(:,WW,cYMO10,t10);
Var_lnWage_WM10=Var_lnWage(:,WM,cYMO10,t10);
dVar_lnWage10=Var_lnWage_WW10-Var_lnWage_WM10;

change_dvar=dVar_lnWage10-dVar_lnWage60;
changelogp=log(relpWW10)-log(relpWW60);

figure(2); figsetup;
plotnamesym(changelogp,change_dvar,ShortNames,9,[],.06,.00);
chadfig2('Change in log p(ww)/p(wm), 1960-2010','Change in Differential Variance(Wage), 1960-2010',1,0);
makefigwide;
%xt=[.01 .04 .16 .5 1 4 16 64];
%xt=[1/64 1/16 1/4 1 4 16 64];
%set(gca,'XTick',log(xt));
%set(gca,'XTickLabel',strmat('.01 .04 .16 .5 1 4 16 64'));
%set(gca,'XTickLabel',strmat('1/64#1/16# 1/4#  1 #  4 # 16 # 64 ','#'));
print VarianceWage_RelPChange.eps
cshow(ShortNames,[changelogp change_dvar],'%8.4f','DlnRelP Change_dVar');
ols(change_dvar,[ones(Noccs,1) changelogp],'ols: Changes in Diff Variance of Wage','Chang_dVar(Wage)',strmat('Constant DlogRelProp'));





