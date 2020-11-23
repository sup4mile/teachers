% SolveEqmBasic_Display.m
%
%   Display equilibrium results -- various statistics including income for Old, Middle, Young by year
%   Called from SolveEqmBasic.m

HYoungModel=squeeze(sum(HModelAll(:,:,1,:),2)); % Noccs x G x YMO x T, like pModel

% Now show the results
fmt='%8.0f %8.0f %8.0f %8.0f %12.3f %8.3f %8.3f %8.3f %12.3f %8.3f %8.3f %8.3f';
tle='Data60 Model60 Data2010 Model Data60 Model Data2010 Model Data60 Model Data2010 Model';
disp ' '; disp ' ';
disp '                                         Wages                          H(i,t)                 HYoung';
youngstuff=[HYoungData(:,1) HYoungModel(:,1) HYoungData(:,6) HYoungModel(:,6)];
stuff=[wData(:,1) wModel(:,1) wData(:,6) wModel(:,6) HData(:,1) HModel(:,1) HData(:,6) HModel(:,6) youngstuff];
cshow(ShortNames,stuff,fmt,tle);


% WM Results
disp ' '; disp ' '; disp '1960 and 2010 quantities';
disp 'White Men p(:,WM) for        Young Cohort                          Middle-Aged Cohort                      Old-Aged Cohort';
tle2='Data60 Model60 Data2010 Model2010 Data60 Model60 Data2010 Model2010 Data60 Model60 Data2010 Model2010';
fmt2='%12.4f %9.4f %9.4f %9.4f %12.4f %9.4f %9.4f %9.4f %12.4f %9.4f %9.4f %9.4f';
%tt=(1:Nyears); cM=7-tt+1;
y1960=1; y2010=6;
cY60=CohortConcordance(y1960,2); cY10=CohortConcordance(y2010,2); Young=1;
cM60=CohortConcordance(y1960,3); cM10=CohortConcordance(y2010,3); Middle=2;
cO60=CohortConcordance(y1960,4); cO10=CohortConcordance(y2010,4); Old=3; % Note: pModel is NxGxYMOxT -- just Y/M/O in the 3rd dimension.
stuffY=[squeeze(pData(:,WM,cY60,y1960)) squeeze(pModel(:,WM,Young,y1960)) squeeze(pData(:,WM,cY10,y2010)) squeeze(pModel(:,WM,Young,y2010))];
stuffM=[squeeze(pData(:,WM,cM60,y1960)) squeeze(pModel(:,WM,Middle,y1960)) squeeze(pData(:,WM,cM10,y2010)) squeeze(pModel(:,WM,Middle,y2010))];
stuffO=[squeeze(pData(:,WM,cO60,y1960)) squeeze(pModel(:,WM,Old,y1960)) squeeze(pData(:,WM,cO10,y2010)) squeeze(pModel(:,WM,Old,y2010))];
cshow(ShortNames,[stuffY stuffM stuffO],fmt2,tle2);

% WW Results
disp ' '; disp ' '; disp '1960 and 2010 quantities';
disp 'White women p(:,WW) for      Young Cohort                          Middle-Aged Cohort                      Old-Aged Cohort';
tle2='Data60 Model60 Data2010 Model2010 Data60 Model60 Data2010 Model2010 Data60 Model60 Data2010 Model2010';
fmt2='%12.4f %9.4f %9.4f %9.4f %12.4f %9.4f %9.4f %9.4f %12.4f %9.4f %9.4f %9.4f';
%tt=(1:Nyears); cM=7-tt+1;
y1960=1; y2010=6;
cY60=CohortConcordance(y1960,2); cY10=CohortConcordance(y2010,2); Young=1;
cM60=CohortConcordance(y1960,3); cM10=CohortConcordance(y2010,3); Middle=2;
cO60=CohortConcordance(y1960,4); cO10=CohortConcordance(y2010,4); Old=3; % Note: pModel is NxGxYMOxT -- just Y/M/O in the 3rd dimension.
stuffY=[squeeze(pData(:,WW,cY60,y1960)) squeeze(pModel(:,WW,Young,y1960)) squeeze(pData(:,WW,cY10,y2010)) squeeze(pModel(:,WW,Young,y2010))];
stuffM=[squeeze(pData(:,WW,cM60,y1960)) squeeze(pModel(:,WW,Middle,y1960)) squeeze(pData(:,WW,cM10,y2010)) squeeze(pModel(:,WW,Middle,y2010))];
stuffO=[squeeze(pData(:,WW,cO60,y1960)) squeeze(pModel(:,WW,Old,y1960)) squeeze(pData(:,WW,cO10,y2010)) squeeze(pModel(:,WW,Old,y2010))];
cshow(ShortNames,[stuffY stuffM stuffO],fmt2,tle2);



% WW Results
disp ' '; disp ' '; disp '1980 and 2010 quantities (the 1960 ones fit by construction)';
disp 'White Women p(:,WW) for      Young Cohort                          Middle-Aged Cohort                      Old-Aged Cohort';
tle2='Data80 Model80 Data2010 Model2010 Data80 Model80 Data2010 Model2010 Data80 Model80 Data2010 Model2010';
fmt2='%12.4f %9.4f %9.4f %9.4f %12.4f %9.4f %9.4f %9.4f %12.4f %9.4f %9.4f %9.4f';
%tt=(1:Nyears); cM=7-tt+1;
y1980=3; y2010=6;
cY80=CohortConcordance(y1980,2); cY10=CohortConcordance(y2010,2); Young=1;
cM80=CohortConcordance(y1980,3); cM10=CohortConcordance(y2010,3); Middle=2;
cO80=CohortConcordance(y1980,4); cO10=CohortConcordance(y2010,4); Old=3; % Note: pModel is NxGxYMOxT -- just Y/M/O in the 3rd dimension.
stuffY=[squeeze(pData(:,WW,cY80,y1980)) squeeze(pModel(:,WW,Young,y1980)) squeeze(pData(:,WW,cY10,y2010)) squeeze(pModel(:,WW,Young,y2010))];
stuffM=[squeeze(pData(:,WW,cM80,y1980)) squeeze(pModel(:,WW,Middle,y1980)) squeeze(pData(:,WW,cM10,y2010)) squeeze(pModel(:,WW,Middle,y2010))];
stuffO=[squeeze(pData(:,WW,cO80,y1980)) squeeze(pModel(:,WW,Old,y1980)) squeeze(pData(:,WW,cO10,y2010)) squeeze(pModel(:,WW,Old,y2010))];
cshow(ShortNames,[stuffY stuffM stuffO],fmt2,tle2);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRODUCTION 
% Look at w*H in data and model 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp ' '; disp ' ';
disp '====================================================================== ';
disp '    PRODUCTION: w*H in data and model for various groups/cohorts';
disp '====================================================================== ';
disp ' ';
TotalLaborProduction=zeros(Ngroups,Ncohorts,Nyears);
YearstoShow=(1:6); %[1 3 6]; % 1960 1980 2010
ShowDetails=0;

for g=1:Ngroups;
    if ShowDetails;
        disp ' ';
        disp '----------------------';
        disp(['   ' GroupNames{g}]);
        disp '----------------------';    
    end;
    for yrs=1:Nyears;
        ProductionYdata =wData(:,yrs) .*squeeze(HAllData(:,g,7-yrs+0,yrs));
        ProductionMdata =wData(:,yrs) .*squeeze(HAllData(:,g,7-yrs+1,yrs));
        ProductionOdata =wData(:,yrs) .*squeeze(HAllData(:,g,7-yrs+2,yrs));
        ProductionYmodel=wModel(:,yrs).*squeeze(HModelAll(:,g,1,yrs));
        ProductionMmodel=wModel(:,yrs).*squeeze(HModelAll(:,g,2,yrs));
        ProductionOmodel=wModel(:,yrs).*squeeze(HModelAll(:,g,3,yrs));
        ProductionData(g,:,yrs) =sum([ProductionYdata ProductionMdata ProductionOdata]);
        ProductionModel(g,:,yrs)=sum([ProductionYmodel ProductionMmodel ProductionOmodel]);
        if ShowDetails;
            disp ' '; disp(Decades(yrs));
            disp '                      Young               Middle               Old';
            tle='Data Model Data Model Data Model';
            fmt='%8.1f %8.1f %12.1f %8.1f %12.1f %8.1f';
            prodstuff=[ProductionYdata ProductionYmodel ProductionMdata ProductionMmodel ProductionOdata ProductionOmodel];
            cshow(ShortNames,prodstuff,fmt,tle);
            disp '------------------------------------------------------------------------';
            cshow('TOTAL          ',sum(prodstuff),fmt,[],[],1);
            disp ' ';
        end;
    end;
end;

for yrs=YearstoShow;
    disp ' '; disp ' '; fprintf('/// %4.0f ///',Decades(yrs));
    disp '          Young                       Middle                      Old';
    tle='Data Model % Data Model % Data Model %';
    fmt='%8.0f %8.0f %8.1f %12.0f %8.0f %8.1f %12.0f %8.0f %8.1f';
    prodstuffY=[ProductionData(:,1,yrs) ProductionModel(:,1,yrs)  100*ProductionModel(:,1,yrs)./ProductionData(:,1,yrs)];
    prodstuffM=[ProductionData(:,2,yrs) ProductionModel(:,2,yrs)  100*ProductionModel(:,2,yrs)./ProductionData(:,2,yrs)];
    prodstuffO=[ProductionData(:,3,yrs) ProductionModel(:,3,yrs)  100*ProductionModel(:,3,yrs)./ProductionData(:,3,yrs)];
    prodstuff=[prodstuffY prodstuffM prodstuffO];
    cshow(GroupNames,prodstuff,fmt,tle);
    disp ' ';
    fprintf('Total income for the year:    Data = %5.0f  Model = %5.0f\n',[sum(sum(ProductionData(:,:,yrs))) sum(sum(ProductionModel(:,:,yrs)))]);
end;





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EARNINGS (after TauW) -- Market Only
% Look at (1-tauw)*w*H in data and model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp ' '; disp ' ';
disp '==================================================================================== ';
disp '    EARNINGS (Market): (1-tauw)*w*H in data and model for various groups/cohorts';
disp '==================================================================================== ';
disp ' ';
TotalLaborIncome=zeros(Ngroups,Ncohorts,Nyears);
ShowDetails=0;

for g=1:Ngroups;
    if ShowDetails;
        disp ' ';
        disp '----------------------';
        disp(['   ' GroupNames{g}]);
        disp '----------------------';    
    end;
    for yrs=1:Nyears;
        IncomeYdata =squeeze(1-TauW(Mkt,g,yrs)).*wData(Mkt,yrs) .*squeeze(HAllData(Mkt,g,7-yrs+0,yrs));
        IncomeMdata =squeeze(1-TauW(Mkt,g,yrs)).*wData(Mkt,yrs) .*squeeze(HAllData(Mkt,g,7-yrs+1,yrs));
        IncomeOdata =squeeze(1-TauW(Mkt,g,yrs)).*wData(Mkt,yrs) .*squeeze(HAllData(Mkt,g,7-yrs+2,yrs));
        IncomeYmodel=squeeze(1-TauW(Mkt,g,yrs)).*wModel(Mkt,yrs).*squeeze(HModelAll(Mkt,g,1,yrs));
        IncomeMmodel=squeeze(1-TauW(Mkt,g,yrs)).*wModel(Mkt,yrs).*squeeze(HModelAll(Mkt,g,2,yrs));
        IncomeOmodel=squeeze(1-TauW(Mkt,g,yrs)).*wModel(Mkt,yrs).*squeeze(HModelAll(Mkt,g,3,yrs));
        IncomeData(g,:,yrs) =nansum([IncomeYdata IncomeMdata IncomeOdata]);
        IncomeModel(g,:,yrs)=nansum([IncomeYmodel IncomeMmodel IncomeOmodel]);
        if ShowDetails;
            disp ' '; disp(Decades(yrs));
            disp '                      Young               Middle               Old';
            tle='Data Model Data Model Data Model';
            fmt='%8.1f %8.1f %12.1f %8.1f %12.1f %8.1f';
            incstuff=[IncomeYdata IncomeYmodel IncomeMdata IncomeMmodel IncomeOdata IncomeOmodel];
            cshow(ShortNames,incstuff,fmt,tle);
            disp '------------------------------------------------------------------------';
            cshow('TOTAL          ',nansum(incstuff),fmt,[],[],1);
            disp ' ';
        end;
    end;
end;

for yrs=YearstoShow;
    disp ' '; disp ' '; fprintf('/// %4.0f ///',Decades(yrs));
    disp '          Young                       Middle                      Old';
    tle='Data Model % Data Model % Data Model %';
    fmt='%8.0f %8.0f %8.1f %12.0f %8.0f %8.1f %12.0f %8.0f %8.1f';
    incstuffY=[IncomeData(:,1,yrs) IncomeModel(:,1,yrs)  100*IncomeModel(:,1,yrs)./IncomeData(:,1,yrs)];
    incstuffM=[IncomeData(:,2,yrs) IncomeModel(:,2,yrs)  100*IncomeModel(:,2,yrs)./IncomeData(:,2,yrs)];
    incstuffO=[IncomeData(:,3,yrs) IncomeModel(:,3,yrs)  100*IncomeModel(:,3,yrs)./IncomeData(:,3,yrs)];
    incstuff=[incstuffY incstuffM incstuffO];
    cshow(GroupNames,incstuff,fmt,tle);
    disp ' ';
    fprintf('Market income for the year:   Data = %5.0f  Model = %5.0f\n',[sum(nansum(IncomeData(:,:,yrs))) sum(sum(IncomeModel(:,:,yrs)))]);
end;



% GDP and Total Labor Income
rho=1-1/sigma;
gdp1=sum((A.*HModel).^rho).^(1/rho);
wHincome=sum(wModel.*HModel);
if any(abs(gdp1-wHincome)>0.1); disp 'Error! GDP is not equal to labor income in SolveEqmBasic. Stopping...'; 
    keyboard; 
else;
    disp 'Production GDP = Total Labor Income (pre-tax) has been confirmed to the nearest dime per person :-).';
end;


disp ' '; disp ' ';
disp 'GDP per person in data and model';
cshow(' ',[Decades GDP YModel ExitFlag],'%12.0f','Year Data Model ExitFlag');

disp ' '; disp ' ';
disp 'GDP (mkt) per person in data and model';
cshow(' ',[Decades GDP_mkt YMktModel],'%12.0f','Year Data Model');

disp ' '; disp ' ';
disp 'GDP per worker (market) in data and model';
cshow(' ',[Decades GDPwkr YwkrModel],'%12.0f','Year Data Model');

disp ' '; disp ' ';
disp 'LFP in data and model';
cshow(' ',[Decades LFPData LFPModel],'%12.0f %12.3f','Year Data Model');

disp ' '; disp ' ';
disp 'Market Wage Gap in data and model';
stuff=[WageGapData(2,:)' WageGapModel(2,:)' WageGapData(3,:)' WageGapModel(3,:)' WageGapData(4,:)' WageGapModel(4,:)'];
cshow(' ',[Decades stuff],'%6.0f %12.3f %8.3f %12.3f %8.3f %12.3f %8.3f','Year WWdata WWmodel BMdata BMmodel BWdata BWmodel');

% disp ' '; disp ' ';
% disp 'Wage Gap in model if everyone worked (Pwork=1)';
% cshow(' ',[Decades WageGapAllModel(2:4,:)'],'%6.0f %12.3f','Year WW BM BW');


disp ' '; disp ' ';
EarningsData= nansum(LaborIncomeReceived_mkt)';
disp 'Earnings (market, after tauw) in data and model';
cshow(' ',[Decades EarningsData EarningsModel],'%12.0f','Year Data Model');

disp ' '; disp ' ';
EarningsData_g=squeeze(nansum(IncomeData,2)); % Sum over cohorts
disp 'Earnings (market, after tauw) by Group in data and model';
stuff=[EarningsData_g(1,:)' EarningsModel_g(1,:)' EarningsData_g(2,:)' EarningsModel_g(2,:)' EarningsData_g(3,:)' EarningsModel_g(3,:)' EarningsData_g(4,:)' EarningsModel_g(4,:)'];
cshow(' ',[Decades stuff],'%6.0f %12.0f %8.0f %12.0f %8.0f %12.0f %8.0f %12.0f %8.0f','Year WMdata WMmodel WWdata WWmodel BMdata BMmodel BWdata BWmodel');


% disp ' '; disp ' ';
% disp 'EarningsAll, model (if everyone worked, i.e. Pwork=1)';
% cshow(' ',[Decades EarningsAllModel],'%12.0f','Year Model');

disp ' '; disp ' ';
disp 'ConsumptionYoung (total) per young person in data and model (and Utility)';
cshow(' ',[Decades ConsumpYoungData ConsumpYoungModel Utility],'%12.0f','Year Data Model Utility');

disp ' '; disp ' ';
disp 'EarningsYoung (market) per young person in data and model';
cshow(' ',[Decades EarningsYoungData EarningsYoungModel],'%12.0f','Year Data Model');

%GDPYoung_Model=sum(wModel.*HYoungModel)';
disp ' '; disp ' '; disp 'GDP (mkt) of only Young Cohorts...';
cshow(' ',[Decades GDPYoung_mkt GDPYoungModel],'%6.0f %8.0f','Year Data Model'); % PerYoung');

% % Home and Total
% disp ' '; disp ' '; disp 'Home and Home+Mkt Production of only Young Cohorts...';
% disp 'Note: HomeY below is per young person, not per young person staying home.';
% HomeShare=HomeEarningsYModel./Home_and_MktOutputY*100;
% cshow(' ',[Decades HomeEarningsYModel Home_and_MktOutputY HomeShare Home_and_MktConsumpY FullConsumpY log(Utility) log(Utility2) log(Util1ms) log(Utilz) log(UtilC) log(UtilC2)],'%6.0f %12.0f %12.0f %12.2f %12.0f %12.0f %12.3f %12.4f %12.4f %12.3f','Year HomeY Home+MktY HomeShare Home_MktC FullCons logUtilityY logUtility2 logU1ms logUtilz logUtilC logUtilC2');

Data =[GDP GDP_mkt GDPwkr LFPData EarningsData ConsumpYoungData EarningsYoungData GDPYoung_mkt WageGapData(2:4,:)' EarningsData_g'];
Model=[YModel YMktModel YwkrModel LFPModel EarningsModel ConsumpYoungModel EarningsYoungModel GDPYoungModel WageGapModel(2:4,:)' EarningsModel_g'];
DataGrowth=log(Data(1,:)./Data(Nyears,:));
ModelGrowth=log(Model(1,:)./Model(Nyears,:));
GrowthShare=(ModelGrowth./DataGrowth)'*100;
tle={'GDP per person','GDP per person (mkt)','GDP (mkt) per worker','Labor Force Participation (LFP)','Earnings','ConsumpYoung (total)','EarningsYoung (market)','GDPYoung (market)','WageGapWW','WageGapBM','WageGapBW','EarningsWM','EarningsWW','EarningsBM','EarningsBW'};
disp ' '; disp ' ';
disp 'Share of growth accounted for by the Full Model:';
cshow(tle,GrowthShare,'%12.1f','Share');

