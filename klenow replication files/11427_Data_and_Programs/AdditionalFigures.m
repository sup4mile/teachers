% AdditionalFigures.m
%
%  Covariance between TauH and TauW for Chang
%  Benchmark levels from Chaining
%  Blau-Kahn female LS elasticities 
%  Charles-Guryan discrimination


clear;
diarychad('AdditionalFigures');
definecolors;

load HowMuchPoorer_Benchmark;

% Covariance between TauH and TauW for Chang
clear stuff;
disp 'Statistics for TW := 1/(1-TauW) and TH := (1+TauH)^eta -- Unweighted'
for g=2:Ngroups;
    disp ' ';
    disp(GroupNames{g});
    for t=1:Nyears;
        tw=1./(1-TauW(2:Noccs,g,t));
        th=(1+TauH(2:Noccs,g,7-t)).^eta;
        CC=corrcoef([tw th]);
        CV=cov([tw th]);
        stuff(t,:)=[Decades(t) CV(1,1) CV(2,2) CV(1,2) CC(1,2)];
    end;
    cshow(' ',stuff,'%8.0f %12.3f','Year Var(TW) Var(TH) Cov(TW,TH) Corr(TW,TH)'); 
end;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDUCATION by year and group.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load SolveEqmBasic_Benchmark;
pModel_Benchmark=pModel; % occ x groups x YMO x year
educ_model=s*25; %  s=1 corresponds to 25 years
for t=1:Nyears;
    ct=7-t;
    EducationModelY(:,t)  =sum(mult(pModel_Benchmark(:,:,1,t),educ_model(:,t))); % YMO for 3rd dim of pModel
    EducationTWTH1960(:,t)=sum(mult(pModel_TWTH(:,:,1,t),educ_model(:,t)));
    EducationDataY(:,t)   =nansum(pData(:,:,ct,t).*Education(:,:,ct,t));
end;

% Note: EducationTWTH1960 is the model solution holding TauW and TauH at 1960 values
EdContributeTWTH=EducationModelY-EducationTWTH1960; % Contribution of TW and TH
EdContributeTWTH=EdContributeTWTH(:,6); % Keep only total contribution 1960-2010

ActualChange=EducationDataY(:,6)-EducationDataY(:,1);
ModelChange =EducationModelY(:,6)-EducationModelY(:,1);
ActualChange_vWM = ActualChange-ActualChange(1);
ModelChange_vWM = ModelChange-ModelChange(1);
EdContributeTWTH_vWM=EdContributeTWTH-EdContributeTWTH(1);

disp ' '; disp ' '; disp ' ';
disp 'Education Results for Young Cohort';
disp '  ModelChg = change in model solution due to all forces';
disp '     TWTH  = change in model solution due to TauW and TauH';
disp ' ';
fprintf('                                                            ---  Relative to WM  ---');
edtle='Data1960 Data2010 Change ModelChg TWTH Change ModelChg TWTH';
cshow(GroupNames,[EducationDataY(:,[1 6]) ActualChange ModelChange EdContributeTWTH ActualChange_vWM ModelChange_vWM EdContributeTWTH_vWM],'%9.2f',edtle);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GDP per person levels graphs
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Benchmark levels from Chaining
%  3/8/19 -- Use Market GDP instead of Home+Mkt = column 2
B=2;
disp ' '; disp ' ';
disp 'Benchmark levels from HowMuchPoorer';
Y_GDP=100*GDP/GDP(1);
Y_Base=100*YBaseline(:,B)/YBaseline(1,B);
Y_twth=100*Y_TWTH(:,B)/Y_TWTH(1,B); % alternative that holds TW/TH constant at 1960 levels

cshow(' ',[Decades Y_Base Y_twth],'%6.0f %8.1f','Year Baseline Y_alt(TWTH)');

figure(1); figsetup; hold on;
%plot(Decades,Y_GDP,'-','Color',myblue);
plot(Decades,Y_Base,'-','Color',myblue,'LineWidth',LW);
%plot(Decades,Y_BothZ,'-','Color',myred);
%plot(Decades,Y_All4,'-','Color',mypurp);
plot(Decades,Y_twth,'-','Color',mygreen,'LineWidth',LW);
vals=(1960:10:2010)';
relabelaxis(vals, num2str(vals),'x');
chadfig2('Year','GDP per person (1960=100)',1,0);
makefigwide;

%text(1985,190,'Data');
%text(1996,190,'Model');
text(1990,210,'Overall');
%text(1998,145,'\tau^h, \tau^w, \it{z}, \Omega_g^{ home}');
%text(2004,118,'\tau^h and \tau^w');
text(1995,150,'\tau^h and \tau^w at 1960 levels');
print('-depsc','BenchmarkLevels');

% plot(Decades,Y_TW,'-','Color',myred,'LineWidth',LW);
% plot(Decades,Y_TH,'-','Color',mypurp,'LineWidth',LW);
% text(2010.5,119,'\tau^h');
% text(2010.5,109,'\tau^w');
% print('-depsc','BenchmarkLevels2');


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % WAGE GAPS levels graphs  9=WW, 10=BM, 11=BW for Wage Gap
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp 'Check to be sure that Wage Gaps are 9=WW, 10=BM, 11=BW in Y_alt=Y_TWTH';
% wait;
% load SolveEqmBasic_Benchmark; % For WageGapBaseline Ngroups x Nyears


% %  Wage Gap levels from Chaining
% FS2=20;
% for g=2:4;
    
%     disp ' '; disp ' ';
%     disp(['Wage Gaps (log) from HowMuchPoorer (constant TWTH=TWTH1960) for ' GroupNames{g}]);
%     WageGap_Data=-100*log(WageGapData(g,:))';
%     WageGap_Model=-100*log(Y_TWTH(:,7+g));
%     cshow(' ',[Decades WageGap_Data WageGap_Model],'%6.0f %8.1f','Year Data Model');

%     figure(1); figsetup; hold on;
%     plot(Decades,WageGap_Data,'-','Color',myblue,'LineWidth',LW);
%     plot(Decades,WageGap_Model,'-','Color',mygreen,'LineWidth',LW);
%     vals=(1960:10:2010)';
%     relabelaxis(vals, num2str(vals),'x');
%     chadfig2('  ','Log points (x100)',1,0,FS2*.7);
%     set(gca,'FontSize',FS2);
%     makefigwide;

%     text(1985,WageGap_Data(3),'Data','FontSize',FS2);
%     text(1990,WageGap_Model(3),'Model','FontSize',FS2); %'\tau^h, \tau^w, \it{z}, \Omega_g^{ home}');
%     wait
%     print('-depsc',['WageGapLevels_' GroupCodes{g}]);
% end;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLAU and KAHN stuff
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Year	Cohort		Blau-Khan	Our Estimate
% 1980	y		0.76	3.13
% 	m		0.75	2.94
% 	o		1.06	3.58
% 1990	y		0.59	1.83
% 	m		0.61	1.74
% 	o		0.68	2.05
% 2000	y		0.33	1.75
% 	m		0.39	1.77
% 	o		0.45	1.64
% Labor supply elasticities for White Women by Year and Cohort
 
% Our results, from CleanandShowTauAZ_Benchmark.log    
%  Year   Young  Middle     Old
%------------------------------

data_hhjk=[  % Updated 1/24/2018
    %  1960    1.48    1.34    1.25
    %1970    1.31    1.21    1.13
  1980    0.92    0.95    1.02
  1990    0.69    0.67    0.75
  2000    0.67    0.68    0.64
    %2010    0.62    0.65    0.63
    ];

data_bk=[
    1980 0.76  0.75  1.06
    1990 0.59  0.61  0.68
    2000 0.33  0.39  0.45
    ];

LSWW_hhjk=vector(data_hhjk(:,2:4));
LSWW_bk  =vector(data_bk(:,2:4));
names=strmat('1980y 1990y 2000y 1980m 1990m 2000m 1980o 1990o 2000o');

figure(1); figsetup;
plotnamesym2(LSWW_bk,LSWW_hhjk,names,10,[],.05,.1); %,1,1);
addolsline(LSWW_bk,LSWW_hhjk,myred,1);
chadfig2('Blau and Kahn elasticities','Model estimates',1,0);
print('-depsc','BlauKahn');
print('-dpsc','AdditionalFigures.ps');



%	Our Measure	Charles Guryan Measure	Weighting Variable
%statefip	weighted_mean_tau	CG_marginal_dis_measure	sample_size_black
data=[
1	1.927	-0.1701046	13701.91
2	1.344	-0.4979347	219.62
4	1.596	-0.6540824	1555.22
5	1.689	-0.4308558	4590.12
6	1.401	-0.6426146	35276.23
8	1.410	-0.6540824	2163.88
9	1.543	-0.5990469	4045.72
10	1.599	-0.4979347	1647.54
12	1.829	-0.5199246	25014.14
13	1.757	-0.363777	25406.71
15	NaN	NaN	211.26
16	NaN	NaN	49.66
17	1.622	-0.5661492	25259.96
18	1.326	-0.6232504	6332.98
19	1.457	-0.8139812	699.75
20	1.528	-0.6412269	2015.34
21	1.478	-0.6426146	3555.34
22	2.000	-0.3841126	17338.21
23	NaN	NaN	95.42
24	1.444	-0.5320681	19695.87
25	1.393	-0.6540824	4643.59
26	1.425	-0.5990469	17583.39
27	1.629	-0.794554	1452.42
28	2.178	-0.363777	11596.53
29	1.501	-0.6362928	7571.32
30	NaN	-0.8139812	33.17
31	2.107	NaN	761.15
32	1.533	NaN	1212.68
33	NaN	-0.66246	136.44
34	1.640	-0.5990469	16638.92
35	NaN	NaN	445.97
36	1.566	-0.5606396	44371.11
37	1.739	-0.3604226	21622.21
38	NaN	-0.66246	37.57
39	1.467	-0.6232504	16808.67
40	1.559	-0.5990469	3075.39
41	1.508	-0.8139812	763.86
42	1.499	-0.6232504	16574.65
44	1.578	-0.6540824	544.74
45	2.104	-0.363777	14584.93
46	NaN	-0.4979347	30.99
47	1.563	-0.5320681	10966.31
48	1.614	-0.5755357	30793.94
49	NaN	-0.66246	188.09
50	NaN	-0.5326762	35.5
51	1.763	-0.5320681	17680.62
53	1.345	-0.6426146	2248.91
54	1.333	-0.8139812	768.18
55	1.681	-0.8139812	3235.31
56	NaN	-0.66246	47.9
];
statenames=strmat('AL AK AZ AR CA CO CT DE FL GA HI ID IL IN IA KS KY LA ME MD MA MI MN MS MO MT NE NV NH NJ NM NY NC ND OH OK OR PA RI SC SD TN TX UT VT VA WA WV WI WY');
statefip=data(:,1);
hhjk=data(:,2);
cg=data(:,3);
weights=data(:,4);


msize=sqrt(weights);
msize=msize/min(msize)*5;

figure(2); figsetup;
s=scatter(hhjk,cg,msize);
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = [0 0.5 0.5];
addolsline(hhjk,cg,myred,1,weights); 
addnames(hhjk,cg,statenames,10,[],.02,.02);
chadfig2('Composite barrier (pooled 1980/90)','Marginal discrimination measure',1,1)
print('-depsc','CharlesGuryan');
print('-dpsc','-append','AdditionalFigures.ps');

diary off;