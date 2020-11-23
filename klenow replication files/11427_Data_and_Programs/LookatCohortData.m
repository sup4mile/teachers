% LookatCohortData.m
%
% Called from ReadCohortData
%  Makes some plots to check that things look okay

% 1990 bw cohort 5 architects missing, whereas small #s in 1970 and 1980!

lname='LookatCohortData';
pname1=[lname '1_' CaseName '.ps'];
pname2=[lname '2_' CaseName '.ps'];
pname3=[lname '3_' CaseName '.ps'];
if exist(pname1); delete(pname1); end;
if exist(pname2); delete(pname2); end;
if exist(pname3); delete(pname3); end;

Occs=[1 8 12 16 23];

% CohortConcordance: c=CohortConcordance(1,2) returns Young in Year1 = cohort 6
CohortConcordance=flipud([
2010 1 2 3
2000 2 3 4
1990 3 4 5
1980 4 5 6 
1970 5 6 7 
1960 6 7 8
]);

% YearConcordance: yr=YearConcordance(6,2) 
%    returns Year in which Cohort 6 is young = Year 1 (1960)
YearConcordance=[
2010 6 % 1 
2000 5 % 2 
1990 4 % 3 
1980 3 % 4  
1970 2 % 5  
1960 1 % 6 
];

clrs='kbrgkbrg';
syms='oxoxoxox';

definecolors;
clrs2=[myred; myblue; mygreen; mypurp; myred; myblue; mygreen; mypurp];

% Aggregate Earnings (for weighting)
AggIncomeData=NumPeople.*Earnings;
AggIncomeData_igt=squeeze(nansum(AggIncomeData,3));
AggIncomeData_it=squeeze(nansum(AggIncomeData_igt,2));
earningsweights=div(AggIncomeData_it,nansum(AggIncomeData_it,1));
earningsweights_avg=mean(earningsweights')';
AggIncomeData_gt=squeeze(nansum(AggIncomeData_igt,1));
earningsweights_gt=div(AggIncomeData_gt,nansum(AggIncomeData_gt,1));
earningsweights_g=mean(earningsweights_gt')';

% EarningsWeights_igt = weights that sum to one, across occupations, for each group in each year
for g=1:Ngroups;
    EarningsWeights_igt(:,g,:)=div(squeeze(AggIncomeData_igt(:,g,:)),squeeze(nansum(AggIncomeData_igt(:,g,:)))');
end;

AggIncomeDataYoung=zeros(Nyears,1);
NumYoung=zeros(Nyears,1);
NumPeoplet=squeeze(nansum(nansum(nansum(NumPeople,3))));
NumPeople_igt=squeeze(nansum(NumPeople,3));
NumPeople_gt=squeeze(nansum(NumPeople_igt,1));

% EARNINGS 
AggEarningsperPerson=squeeze(nansum(nansum(nansum(AggIncomeData))))./NumPeoplet;
disp 'Total earnings per person';
cshow(' ',[Decades AggEarningsperPerson],'%8.0f %12.0f','Decade Earnings');
Earn=AggEarningsperPerson;

disp ' ';
disp 'Average growth rates of Earnings';
disp '---------------------------';
g=[log(Earn(2:end)./Earn(1:(end-1)))]./(Decades(2:end)-Decades(1:(end-1)));
gAll=log(Earn(end)/Earn(1))/(2010-1960);
yr0=Decades; yr0(end)=1960; yrT=[Decades(2:end); 2010];
cshow(' ',[yr0 yrT [g; gAll]],'%6.0f %5.0f %9.4f');

% EARNINGS by group
AggEarningsperPerson_g=(AggIncomeData_gt./NumPeople_gt)';
disp ' ';
disp 'Total earnings per person by group';
cshow(' ',[Decades AggEarningsperPerson_g],'%8.0f %12.0f','Decade WM WW BM BW');
Earn_g=AggEarningsperPerson_g;

disp ' ';
disp 'Average growth rates of Earnings by Group';
disp '---------------------------';
g=[log(Earn_g(2:end,:)./Earn_g(1:(end-1),:))]./(Decades(2:end)-Decades(1:(end-1)));
gAll_g=log(Earn_g(end,:)./Earn_g(1,:))/(2010-1960);
yr0=Decades; yr0(end)=1960; yrT=[Decades(2:end); 2010];
tle='Year0 YearT WM WW BM BW';
cshow(' ',[yr0 yrT [g; gAll_g]],'%6.0f %5.0f %9.4f %9.4f %9.4f %9.4f',tle);

disp ' ';
disp '================================';
disp 'BACK OF THE ENVELOPE CALCULATION';
disp '================================';
disp ' ';
fprintf('Growth rate of all people: %8.4f\n',gAll);
fprintf('Growth rate of    WM     : %8.4f\n',gAll_g(1));
BOFEshare=1-gAll_g(1)/gAll;
fprintf('Rough share of growth from closing gaps: %8.3f\n',BOFEshare);
disp ' ';


% EARNINGS YOUNG
for t=1:Nyears;
    AggIncomeDataYoung(t)=nansum(nansum(AggIncomeData(:,:,7-t,t)));
    NumYoung(t)=nansum(nansum(NumPeople(:,:,7-t,t)));
end;
disp ' ';
disp 'Earnings of Young: Per total population and Per person';
cshow(' ',[Decades AggIncomeDataYoung./NumPeoplet AggIncomeDataYoung./NumYoung NumYoung./NumPeoplet*100],'%8.0f %12.0f %12.0f %12.1f','Year PerTotalPop PerYoung Young/Pop');

for i=1:length(Occs);
    for g=1:4; % Groups
        tle=[ShortNames{Occs(i)} ': ' GroupNames{g}];

        % Now plot all cohorts over time for iSec
        sfigure(1); figsetup; hold on;
        for c=1:8;
            % Find the years for a given cohort
            yrs=Decades(any(CohortConcordance'==c)); 
            [tf yy]=ismember(yrs,Decades);
            plot(yrs,squeeze(p(Occs(i),g,c,yy)),syms(c),'Color',clrs2(g,:)); %,cat(2,clrs(c),syms(c)));
            plot(yrs,squeeze(p(Occs(i),g,c,yy)),cat(2,clrs(c),'-'));
        end;
        title(tle);
        chadfig2('  ','p(.)',1,0);
        makefigwide;
        print('-dpsc','-append',pname1);
    end;
end;


% Relative propensity
relp=zeros(size(p));
for g=1:Ngroups;
    relp(:,g,:,:)=p(:,g,:,:)./p(:,1,:,:);
end;



for i=1:length(Occs);
    for g=1:4; % Groups
        tle=[ShortNames{Occs(i)} ': ' GroupNames{g}];

        % Now plot all cohorts over time for iSec
        sfigure(1); figsetup; hold on;
        for c=1:8;
            % Find the years for a given cohort
            yrs=Decades(any(CohortConcordance'==c)); 
            [tf yy]=ismember(yrs,Decades);
            plot(yrs,squeeze(relp(Occs(i),g,c,yy)),cat(2,clrs(c),syms(c)));
            plot(yrs,squeeze(relp(Occs(i),g,c,yy)),cat(2,clrs(c),'-'));
        end;
        chadfig2('  ','Rel. propensity, p(.)/p(wm)',1,0);
        makefigwide;
        title(tle);
        print('-dpsc','-append',pname2);
    end;
end;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WageBar as Geometric Average
%  -- from the Arithmetic average variable "Earnings"
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('Wage');
    Wage_controls=Wage;
    WageNominal_controls=WageNominal;
    clear Wage WageNominal % To avoid confusion; Wage_control is a regression residual
end;
%WageBar=GammaBar/GammaBase*Earnings; % Using *inflated* Arithmetic mean for WageBar for now.

% 2019: Now using Erik's geometric mean for WageBar
WageBar=Earn_geo;

% Wage Gap -- Ratio of WageBar
%   Graph of earnings-weighted average for each group

wagegap=zeros(size(p));
for g=1:Ngroups;
    %    wagegap(:,g,:,:)=Wage(:,g,:,:)./Wage(:,WM,:,:);
    wagegap(:,g,:,:)=WageBar(:,g,:,:)./WageBar(:,WM,:,:);
    for c=1:Ncohorts;
        wagegap_gct(g,c,:)=nansum(mult(squeeze(wagegap(:,g,c,:)),earningsweights_avg));
    end;
end;


% Wage Gap Graph (time vs cohort effects ==> clear cohort effect!)
% (analogous to Erik's graph on slide 49 from 10/26/17
GroupCodes={'WM','WW','BM','BW'};
for g=2:4; % Groups
    tle=GroupNames{g};
    % Now plot all cohorts over time
    sfigure(1); figsetup; hold on;
    for c=1:8;
        % Find the years for a given cohort
        yrs=Decades(any(CohortConcordance'==c)); 
        [tf yy]=ismember(yrs,Decades);
        plot(yrs,log(squeeze(wagegap_gct(g,c,yy))),syms(c),'Color',clrs2(c,:)); %,cat(2,clrs(c),syms(c)));
        plot(yrs,log(squeeze(wagegap_gct(g,c,yy))),'-','Color',clrs2(c,:)); 
    end;
    %ax=axis; ax(4)=0; axis(ax);
    %if g==2;
    %    vals=(-1:.2:0)';
    %    relabelaxis(vals, num2str(vals));
    %end;
    vals=(1960:10:2010)';
    relabelaxis(vals, num2str(vals),'x');
    chadfig2('  ','Log Wage Gap',1,0);
    makefigwide;
    print('-depsc',['WageGap' GroupCodes{g} '_' CaseName]);
    title(tle,'FontName','Helvetica','FontSize',14);
    print('-dpsc','-append',pname2);
end;



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wage Gaps versus Relative Propensity Graphs (as in Version 3.0)
%  -- use young cohort and the Earnings = WageBar series from Erik
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Occupational wage gaps for white women in 1980 vs relative p's
ymo=1+1; % Young
t80=find(Decades==1980);
cYMO80=CohortConcordance(t80,ymo);
relpWW80=relp(:,WW,cYMO80,t80);
%occwagegaps1980=-log(Wage_controls(:,WW,cYMO80,t80)./Wage_controls(:,WM,cYMO80,t80));
occwagegaps1980=-log(WageBar(:,WW,cYMO80,t80)./WageBar(:,WM,cYMO80,t80));

figure(1); figsetup;
plotnamesym2(log(relpWW80),occwagegaps1980,ShortNames,9,[],.06,.00);
chadfig2('Relative propensity, p(ww)/p(wm)','Occupational wage gap (logs)',1,0);
makefigwide;
%xt=[.01 .04 .16 .5 1 4 16 64];
xt=[1/64 1/16 1/4 1 4 16 64];
set(gca,'XTick',log(xt));
%set(gca,'XTickLabel',strmat('.01 .04 .16 .5 1 4 16 64'));
set(gca,'XTickLabel',strmat('1/64#1/16# 1/4#  1 #  4 # 16 # 64 ','#'));
if HighQualityFigures; wait; end;
print('-depsc',['OccupationalGaps1980a_' CaseName]);
%print OccupationalGaps1980a.eps
cshow(ShortNames,[relpWW80 occwagegaps1980],'%8.4f','RelP Gap');
ols(occwagegaps1980,[ones(Noccs,1) log(relpWW80)],'ols: Levels','logwagegap',strmat('Constant logRelProp'));

% Now the change for Figure 2
t60=find(Decades==1960);
cYMO60=CohortConcordance(t60,ymo);
relpWW60=relp(:,WW,cYMO60,t60);
%occwagegaps1960=-log(Wage_controls(:,WW,cYMO60,t60)./Wage_controls(:,WM,cYMO60,t60));
occwagegaps1960=-log(WageBar(:,WW,cYMO60,t60)./WageBar(:,WM,cYMO60,t60));

t10=find(Decades==2010);
cYMO10=CohortConcordance(t10,ymo);
relpWW10=relp(:,WW,cYMO10,t10);
%occwagegaps2010=-log(Wage_controls(:,WW,cYMO10,t10)./Wage_controls(:,WM,cYMO10,t10));
occwagegaps2010=-log(WageBar(:,WW,cYMO10,t10)./WageBar(:,WM,cYMO10,t10));

changewage=occwagegaps2010-occwagegaps1960;
changelogp=log(relpWW10)-log(relpWW60);

figure(2); figsetup;
plotnamesym2(changelogp,changewage,ShortNames,9,[],.06,.00);
chadfig2('Change in log p(ww)/p(wm), 1960-2010','Change in log wage gap, 1960-2010',1,0);
makefigwide;
%xt=[.01 .04 .16 .5 1 4 16 64];
%xt=[1/64 1/16 1/4 1 4 16 64];
%set(gca,'XTick',log(xt));
%set(gca,'XTickLabel',strmat('.01 .04 .16 .5 1 4 16 64'));
%set(gca,'XTickLabel',strmat('1/64#1/16# 1/4#  1 #  4 # 16 # 64 ','#'));
if HighQualityFigures; wait; end;
print('-depsc',['OccupationalGaps1980b_' CaseName]);
%print OccupationalGaps1980b.eps
cshow(ShortNames,[changelogp changewage],'%8.4f','DlnRelP DlnGap');
ols(changewage,[ones(Noccs,1) changelogp],'ols: Changes','Dlogwagegap',strmat('Constant DlogRelProp'));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Var_lnIncome and Var_lnWage versus Relative Propensity Graphs 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isequal(CaseName,'Benchmark');
    Variance_IncomeWages
end;


% tauhat -- wage gap by Occupation/group/cohort/year...
%   
%   relp = tauhat^(-theta) * [wage(g)/wage(wm)]^(-theta(1-eta))
%
%   tauhat = relp^(-1/theta) * [wage(g)/wage(wm)]^(-(1-eta))



% tauhat=zeros(size(p));
% for i=1:Noccs;
%     tauhat(i,:,:,:)=squeeze(relp(i,:,:,:)).^(-1/theta) .* wagegap.^(-(1-eta));
% end;

disp ' '; disp ' ';
disp 'Computing tauhat for *all* c,t:   tauhat=relp.^(-(1-dlta)/theta) .* wagegap.^(-(1-eta));'
disp 'just to see what they look like. We only use the Young (c,c) version.';
disp 'This version ignores the Tbar/gamma stuff. But for baseline case and Young = correct.';
disp ' ';
tauhat_all=relp.^(-(1-dlta)/theta) .* wagegap.^(-(1-eta));
tauhat_all(1,:,:,:)=1; % Fix HOME

for i=1:length(Occs);
    for g=1:4; % Groups
        tle=[ShortNames{Occs(i)} ': ' GroupNames{g}];

        % Now plot all cohorts over time for iSec
        sfigure(1); figsetup; hold on;
        for c=1:8;
            % Find the years for a given cohort
            yrs=Decades(any(CohortConcordance'==c)); 
            [tf yy]=ismember(yrs,Decades);
            plot(yrs,squeeze(tauhat_all(Occs(i),g,c,yy)),cat(2,clrs(c),syms(c)));
            plot(yrs,squeeze(tauhat_all(Occs(i),g,c,yy)),cat(2,clrs(c),'-'));
        end;
        title(tle);
        chadfig2('  ','tauhat',1,0);
        makefigwide;
        print('-dpsc','-append',pname3);
    end;
end;

% Finally, create key graphs for the paper: key occupations, young tauhat for WW and BM
% First, extract the young tauhat
onames={'CompositeTauOccWM','CompositeTauOccWW','CompositeTauOccBM','CompositeTauOccBW'};
if ~exist('WideFigures'); WideFigures=1; end;
if ~exist('HighQualityFigures'); HighQualityFigures=0; end;
for t=1:Nyears;
    ct=7-t;
    tauhat(:,:,t)=squeeze(tauhat_all(:,:,ct,t)); % Noccs x Ngroups x Nyears.  Previously tauhat_y
end;
tauhat_y=tauhat; % legacy

% Copying code from EstimateTauZ to fix missing tauhat_y for use in NoBrawny case
% (where we use tauhat_y to infer T(i,g) to eliminate frictions)

    % Fix missing values in a systematic way
    disp ' '; disp '===========================================================';
    disp 'Fixing tauhat missing values --> tauhat_y_cleaned (just for robustness cases)';
    disp ' -- Use first non-missing value';
    disp ' ';  disp '===========================================================';

    % tauhat_y
    tauhat_y_cleaned=tauhat;
    for g=2:Ngroups; disp ' ';
        disp '*********************';
        disp(GroupNames{g});
        disp '*********************'; disp ' ';
        for i=2:Noccs;
            % tauhat_y
            missingtauhat_y_cleaned=(isnan(squeeze(tauhat_y_cleaned(i,g,:))) | isinf(squeeze(tauhat_y_cleaned(i,g,:))));
            if all(missingtauhat_y_cleaned); 
                disp 'All tauhat_y_cleaned missing. Stopping...'; keyboard;
            elseif any(missingtauhat_y_cleaned);
                cshow(['tauhat_y_cleaned ' ShortNames{i}],squeeze(tauhat_y_cleaned(i,g,:))','%8.4f',[],'nonee',1);
                tauhat_y_cleaned(i,g,:)=fixmissing(squeeze(tauhat_y_cleaned(i,g,:)),missingtauhat_y_cleaned);
                cshow(['tauhat_y_cleaned ' ShortNames{i}],squeeze(tauhat_y_cleaned(i,g,:))','%8.4f',[],'nonee',1);
            end;
        end; % Occupations
    end; % Groups




occExe=1+1;
occDoc=7+1;
occTea=11+1;
occLaw=15+1;
occSal=21+1;
occSec=22+1;
occFoo=36+1;
occCst=47+1;
myoccs=[occDoc occLaw occSec occCst occTea];
definecolors


% Doctors Lawyer Home Secretaries
%ytle='{\bf Barrier measure, $$\hat\tau$$}';
%ytle='Composite barrier, \tau';
ytle='Composite barrier';
if HighQualityFigures; ytle=' '; end;
for g=2:Ngroups;
  sfigure(g); figsetup;
  plot(Decades,log(vector(tauhat(occDoc,g,:))),'Color',myblue,'LineWidth',LW); %
  plot(Decades,log(vector(tauhat(occLaw,g,:))),'Color',mygreen,'LineWidth',LW); %
  plot(Decades,log(vector(tauhat(occSec,g,:))),'Color',mypurp,'LineWidth',LW); %
  plot(Decades,log(vector(tauhat(occCst,g,:))),'Color',myorng,'LineWidth',LW); %
  plot(Decades,log(vector(tauhat(occTea,g,:))),'Color',myred,'LineWidth',LW); %
  relabelaxis(log([1/8 1/4 1/2 1 2 4 8 16 32]),'1/8 1/4 1/2 1 2 4 8 16 32','y');
  disp(GroupNames{g})
  cshow(ShortNames(myoccs),squeeze(tauhat_y(myoccs,g,:)),'%8.3f','1960 1970 1980 1990 2000 2010');

  %  if g==2; 
  %  set(gca,'YTick',[0 1 2 3 4 5]);
  %  set(gca,'YTickLabel',strmat('0 1 2 3 4 5'));
  %end;
  vals=(1960:10:2010)';
  relabelaxis(vals, num2str(vals),'x');
  chadfig2('  ',ytle,1,0);
  if WideFigures; makefigwide; end;
%  title(tle{g});
  if HighQualityFigures; 
    putstr('Doctors'); %,myblue);
    putstr('Lawyers'); %,mygreen);
    putstr('Secretaries'); %,mypurp);
    putstr('Construction'); %,myorng);
    putstr('Teachers'); %,myred);
    wait; 
  end;  % For adjusting years...
  print('-depsc',[onames{g} '_' CaseName]);
end;




    
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot mean and var of composite tau
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
definecolors
meannames='MeanTauHat';
varnames='VarTauHat';
relpnames='StdlnRelP';
gapnames='StdWageGap';
Mkt=2:Noccs;


for t=1:Nyears;
  meantau(t,:)=nansum(mult(tauhat_y_cleaned(:,:,t),earningsweights_avg));
  meanlogtau(t,:)=nansum(mult(log(tauhat_y_cleaned(:,:,t)),earningsweights_avg));
  diff=chadminus(log(tauhat_y_cleaned(:,:,t)),meanlogtau(t,:));
  varlogtau(t,:)=nansum(mult(diff.^2,earningsweights_avg));
  
  % std of log(relp)
  lnrelp_young(:,:,t)=log(relp(:,:,7-t,t));
  lnrelp_young(isinf(lnrelp_young))=NaN;
  lnrelp_young(isnan(lnrelp_young))=log(0.001); %min(nanmin(lnrelp_young(:,:,t)));
  meanrelp=nansum(mult(lnrelp_young(:,:,t),earningsweights_avg));
  diff=chadminus(lnrelp_young(:,:,t),meanrelp);
  std_relp_young(t,:)=sqrt(nansum(mult(diff.^2,earningsweights_avg)));  
  
  % std of gap := log(wage ratio)
  wage_young(:,:,t)=Earnings(:,:,7-t,t);
  gap(:,:,t)=log(div(wage_young(:,:,t),wage_young(:,WM,t)));
  gap(isinf(gap))=NaN;
  missinggap=isnan(gap(:,:,t)); missinggap(1,:)=0;
  meangap=nansum(mult(gap(:,:,t),earningsweights_avg));
  gap(missinggap(:,WW),WW,t)=meangap(:,WW); %nanmean(gap(~missinggap(:,WW),WW,t));
  gap(missinggap(:,BM),BM,t)=meangap(:,BM); %nanmean(gap(~missinggap(:,BM),BM,t));
  gap(missinggap(:,BW),BW,t)=meangap(:,BW); %nanmean(gap(~missinggap(:,BW),BW,t));
  diff=chadminus(gap(:,:,t),meangap);
  std_wagegap(t,:)=sqrt(nansum(mult(diff.^2,earningsweights_avg)));
end;


FS2=20 % Font Size for 2x2 Figures
sfigure(1); figsetup;
if WideFigures; makefigwide; end;
%plot(Decades,meantau(:,1),'b-'); hold on;
plot(Decades,meantau(:,WW),'Color',myblue,'LineWidth',LW); hold on;
plot(Decades,meantau(:,BM),'Color',mygreen,'LineWidth',LW);
plot(Decades,meantau(:,BW),'Color',mypurp,'LineWidth',LW);
vals=(1960:10:2010)';
relabelaxis(vals, num2str(vals),'x');
ytle='Mean (weighted) across occupations';
if HighQualityFigures; ytle=' '; end;
chadfig2('  ',ytle,1,0,FS2*.7);
set(gca,'FontSize',FS2);
if 1; %HighQualityFigures;
  text(1963,4,mlstring('White\\    women'),'FontSize',FS2); 
  text(1972,2.3,'Black men','FontSize',FS2); %,'Color',mygreen);
  text(1971,7,mlstring('Black women'),'FontSize',FS2); 
end;
print('-depsc',[meannames '_' CaseName]);
title('Mean of Composite TauHat');
print('-dpsc','-append',pname3);


sfigure(2); figsetup;
%plot(Decades,meantau(:,1),'b-'); hold on;
plot(Decades,varlogtau(:,WW),'Color',myblue,'LineWidth',LW); hold on;
plot(Decades,varlogtau(:,BM),'Color',mygreen,'LineWidth',LW);
plot(Decades,varlogtau(:,BW),'Color',mypurp,'LineWidth',LW);
vals=(1960:10:2010)';
relabelaxis(vals, num2str(vals),'x');
ytle='Variance (weighted) of log';
if HighQualityFigures; ytle=' '; end;
chadfig2('  ',ytle,1,0,FS2*.7);
set(gca,'FontSize',FS2);
if WideFigures; makefigwide; end;
if 1; % HighQualityFigures;
  text(1977,.55,'White women','FontSize',FS2); %,'Color',myblue);
  text(1980,.2,'Black men','FontSize',FS2); %,'Color',mygreen);
  text(1990,.8,'Black women','FontSize',FS2); %,'Color',mypurp);
end;
print('-depsc',[varnames '_' CaseName]);
title('Variance of log tauhat');
print('-dpsc','-append',pname3);


% Now we make the figures that Erik wants at the start of the paper

% stdev of relative propensities
sfigure(3); figsetup;
greygrid(Decades,(0.5:.5:2));
plot(Decades,std_relp_young(:,WW),'Color',myblue,'LineWidth',LW); hold on;
plot(Decades,std_relp_young(:,BM),'Color',mygreen,'LineWidth',LW); 
plot(Decades,std_relp_young(:,BW),'Color',mypurp,'LineWidth',LW); 
vals=(0:.5:2.5)';
relabelaxis(vals,num2str(vals));
vals=(1960:10:2010)';
relabelaxis(vals, num2str(vals),'x');
chadfig2('  ','  ',1,0);
if WideFigures; makefigwide; end;
if 1; %HighQualityFigures;
  text(1996,1.35,'White women'); %,'Color',myblue);
  text(1980,0.75,'Black men'); %,'Color',mygreen);
  text(1967,2.3,'Black women'); %,'Color',mypurp);
end;
print('-depsc',[relpnames '_' CaseName]);
tle='Stdev of log(relp) across occupations';
title(tle);
print('-dpsc','-append',pname3);
disp ' '; disp(tle);
cshow(' ',[Decades std_relp_young],'%6.0f %8.3f','Year WM WW BM BW');

% stdev of Wage Gaps
sfigure(4); figsetup;
greygrid(Decades,(0.1:.1:.4));
plot(Decades,std_wagegap(:,WW),'Color',myblue,'LineWidth',LW); hold on;
plot(Decades,std_wagegap(:,BM),'Color',mygreen,'LineWidth',LW); 
plot(Decades,std_wagegap(:,BW),'Color',mypurp,'LineWidth',LW); 
vals=(0:.1:.5)';
relabelaxis(vals,num2str(vals));
vals=(1960:10:2010)';
relabelaxis(vals, num2str(vals),'x');
chadfig2('  ','  ',1,0);
if WideFigures; makefigwide; end;
if 1; %HighQualityFigures;
  text(2000,.05,'White women'); %,'Color',myblue);
  text(1963,.12,'Black men'); %,'Color',mygreen);
  text(1975,.26,'Black women'); %,'Color',mypurp);
end;
print('-depsc',[gapnames '_' CaseName]);
tle='Stdev of log(Wage Gap) across occupations';
title(tle);
print('-dpsc','-append',pname3);
disp ' '; disp(tle);
cshow(' ',[Decades std_wagegap],'%6.0f %8.3f','Year WM WW BM BW');

   
    
% Note, clearing stuff...
clear data*
clear c Cohort i decindx fmt fname g group occnum tf tle total xNumPeople year yrs yy
save(['CohortData_' CaseName]);