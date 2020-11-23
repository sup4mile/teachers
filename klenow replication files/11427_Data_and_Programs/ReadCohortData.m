% ReadCohortData
%
% Reads the cohort data assembled by Erik in chad_output_file_[Date].csv
% (see the .xls file for more details).



if ~exist('CaseName');
    CaseName=input('What would you like to call the Case Name? ','s');
end;
diarychad('ReadCohortData',CaseName);
Names67Occupations; % load the names and "brawny" occupation index.
ShowParameters;

%Year,Group,Cohort,Occupation Number,Weighted Total People in Occ,Avg Occ Income,Education,Wage
%2010,10,1,0,7610209.5,,12.7,
%2010,10,1,1,1965792.0,54377.4,15.0,

% %year,group,cohort,occ_code,num,occ_income,occ_grade,occ_wage
% 2010,0,0,0,2.43E+07,3208.442,12.60451,
% 2010,0,0,1,8273344,82963.73,15.10609,
% 2010,0,0,2,4476900,68100.81,15.21636,

% Version 6 code
fname='chad_output_file_2019_01_24.csv'
fmt='%f %f %f %f %f %f %f %f %f %f %f %f %f %f';
[year,group,cohort,occnum,dataNumPeople,dataEarnings,dataEducation,dataWage,dataLnEarn_arith,dataLnEarn_geo,dataLnWage_arith,dataLnWage_geo,datavar_lninc,datavar_lnwage]...
   =textread(fname,fmt,'headerlines',1,'delimiter',',','emptyvalue',NaN);
year(year==2012)=2010; % In latest version of Erik's csv file, he records year as 2012

occnum=occnum+1;  % So Home=1 rather than Home=0

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize key variables
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nrecords=length(dataNumPeople);
Noccs=67
Ngroups=4  % wm ww bm bw
Ncohorts=8
Decades=[1960 1970 1980 1990 2000 2010]';
Nyears=length(Decades);

NumPeople=zeros(Noccs,Ngroups,Ncohorts,Nyears);
Education=NumPeople;
EarningsNominal=NumPeople;
WageNominal=NumPeople;
Earnings=NumPeople;
Wage=NumPeople;
Earn_arith=NumPeople; % Should be same as Earnings (checked)
Earn_geo=NumPeople;   % Geometric average rather than arithmetic average
Wage_arith=NumPeople; % Not using for now (includes controls)
Wage_geo=NumPeople;   % Not using for now (includes controls)
Var_lnIncome=NumPeople;
Var_lnWage=NumPeople;

disp 'Putting earnings / wage in $2009 constant using PCE Deflator';
% See PCEDeflatorNIPA.txt from https://research.stlouisfed.org/fred2/series/DPCERD3A086NBEA
pce=[
   17.535   % 1960-01-01
   22.311   % 1970-01-01
   43.959   % 1980-01-01
   67.440   % 1990-01-01
   83.131   % 2000-01-01
  101.653   % 2010-01-01
]';


% Totals across all groups (provided by Erik in the .xls file)
AllNumPeople=zeros(Noccs,Ncohorts,Nyears);
AllEducation=AllNumPeople;
AllEarnings=AllNumPeople;



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reshape into the multidimensional matrices
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:Nrecords;
    if cohort(i)~=0; % Erik added this for *all* cohorts in a single year (just ignore it)
        if group(i)==0; % All -- new code as of 7/25/16
            AllNumPeople(occnum(i),cohort(i),find(Decades==year(i)))=dataNumPeople(i);
            AllEarnings(occnum(i),cohort(i),find(Decades==year(i)))=dataEarnings(i);
            AllEducation(occnum(i),cohort(i),find(Decades==year(i)))=dataEducation(i);
        else; % One of our 4 key groups
            if year(i)>1950; % Ignore the 1950 data for now...
                decindx=find(Decades==year(i)); % e.g. decindx=2 for year=1970
                NumPeople(occnum(i),group(i),cohort(i),decindx)=dataNumPeople(i);
                if dataEarnings(i)==0; dataEarnings(i)=NaN; end; % If 0 earnings
                EarningsNominal(occnum(i),group(i),cohort(i),decindx)=dataEarnings(i);
                Earnings(occnum(i),group(i),cohort(i),decindx)=dataEarnings(i)/pce(decindx)*pce(Nyears);
                Education(occnum(i),group(i),cohort(i),decindx)=dataEducation(i);
                WageNominal(occnum(i),group(i),cohort(i),decindx)=dataWage(i);
                Wage(occnum(i),group(i),cohort(i),decindx)=dataWage(i)/pce(decindx)*pce(Nyears);
                Earn_arith(occnum(i),group(i),cohort(i),decindx)=exp(dataLnEarn_arith(i))/pce(decindx)*pce(Nyears);
                Earn_geo(occnum(i),group(i),cohort(i),decindx)=exp(dataLnEarn_geo(i))/pce(decindx)*pce(Nyears);
                Var_lnIncome(occnum(i),group(i),cohort(i),decindx)=datavar_lninc(i);
                Var_lnWage(occnum(i),group(i),cohort(i),decindx)=datavar_lnwage(i);
            end;        
        end;
    end; % if cohort(i)~=0
end;

% Shut off any earnings in the Home sector -- remnants of Erik's program 7/25/16
AllEarnings(1,:,:,:)=NaN;
EarningsNominal(1,:,:,:)=NaN;
Earnings(1,:,:,:)=NaN;
Earn_arith(1,:,:,:)=NaN;
Earn_geo(1,:,:,:)=NaN;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute p == fraction of WW in Cohort 3 in 2000 who are lawyers
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p=zeros(size(NumPeople));
% Treat "NaN as 0 for purposes of computing p
xNumPeople=NumPeople;
xNumPeople(isnan(xNumPeople))=0;
total=sum(xNumPeople,1); % Add across occupations
for i=1:Noccs;
    p(i,:,:,:)=xNumPeople(i,:,:,:)./total;
end;


% pDataYWM -- p in the data for Young WM
for t=1:Nyears;
    pDataYWM(:,t)=p(:,1,7-t,t);
end;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute q(g,c,t) == fraction of Population who are WW in Cohort 3 in 2000 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xNumPeople_t=sum(squeeze(sum(sum(xNumPeople)))); % 1xT vector date t "populations"
xNumPeople_gct=squeeze(sum(xNumPeople));         % Sum over occupations
q=zeros(Ngroups,Ncohorts,Nyears)*NaN;
for t=1:Nyears;
    q(:,:,t)=squeeze(xNumPeople_gct(:,:,t))/xNumPeople_t(t);
end;




% Adjust Earnings if WageGapAdjustmentFactor=1/2 or Zero.
%  That is, Earnings = (1-WageGapAdjustmentFactor)*Earnings(WM) + WageGapAdjustmentFactor*Earnings(g)
%   Zero ==> Earnings = WM earnings, so no wage gap.
%   1/2  ==>  equally-weighted average of own and WM earnings.
for g=1:Ngroups;
    Earnings(:,g,:,:)=(1-WageGapAdjustmentFactor)*Earnings(:,WM,:,:) + WageGapAdjustmentFactor*Earnings(:,g,:,:);
    Wage(:,g,:,:)=(1-WageGapAdjustmentFactor)*Wage(:,WM,:,:) + WageGapAdjustmentFactor*Wage(:,g,:,:);
    Earn_arith(:,g,:,:)=(1-WageGapAdjustmentFactor)*Earn_arith(:,WM,:,:) + WageGapAdjustmentFactor*Earn_arith(:,g,:,:);
    Earn_geo(:,g,:,:)=(1-WageGapAdjustmentFactor)*Earn_geo(:,WM,:,:) + WageGapAdjustmentFactor*Earn_geo(:,g,:,:);
end;
Wage(1,:,:,:)=NaN;  % No home data. Erik says to 'NaN' this out bc of "pollution" from his program.


% Education YWM for calibrating phiFARM
EducationYWM=zeros(Noccs,Nyears);
for t=1:Nyears;
    EducationYWM(:,t)=Education(:,1,7-t,t);
end;
disp ' ';
disp 'Education of Young WM';
cshow(ShortNames,EducationYWM,'%8.4f','1960 1970 1980 1990 2000 2010');
disp ' '

% disp ' ';
% fprintf(['OccupationtoIdentifyPhi = %2.0f ' OccupationNames(OccupationtoIdentifyPhi,:) '\n'],OccupationtoIdentifyPhi);

% % Allow use of FARM explicitly
% % (comment out next "if" line when using sales.)
% %if OccupationtoIdentifyPhi~=42;
%     disp 'Implied values of Phi for this occupation:'
%     sKey=EducationYWM(OccupationtoIdentifyPhi,:)/25; % 25 years is denominator
%     PhiKeyOcc=(1-eta)/beta.*sKey./(1-sKey)
%     %else;
%     %PhiKeyOcc=phiFarm
%     %end;

%     %wait;

% Now call LookatCohortData to create some graphs
% and save the data for estimation.




LookatCohortData

diary off;
