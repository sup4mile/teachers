% SolveEqmBasic.m    
%
%  Takes the TalentData.mat data on Taus, Z's, A(i,t).
%  Solves for the equilibrium w(i) and Y
%
%  See 2015-06-02-SolvingGE.pdf notes (old but somewhat useful?)
%
%  Method: For each year,
%    1. Guess values for {mgtilde}, Y ==> 5 unknowns
%    2. Solve for {wi} from Hi^supply = Hi^demand
%    3. Compute mghat, Yhat 
%    4. Iterate until converge.

clear; global CaseName;
diarychad('SolveEqmBasic',CaseName);

global Noccs Ngroups Ncohorts Nyears CohortConcordance TauW_Orig pData HAllData q ShortNames
global TauW_C phi_C mgtilde_C w_C GeoArith_adjust mgEstimated % For keeping track of history in solution

load(['TalentData_' CaseName]); % From EstimateTauZ and earlier programs
ShowParameters;

[YModel,YMktModel,EarningsModel,YwkrModel,LFPModel,ConsumpYoungModel,EarningsYoungModel,GDPYoungModel,WageGapModel,EarningsModel_g,Utility,wModel,HModel,HModelAll,pModel,ExitFlag]=SolveForEqm(TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar);
% Now show results using a separate program (so we can call it elsewhere when testing)
SolveEqmBasic_Display


GDPBaseline=YModel;
GDPMktBaseline=YMktModel;
GDPwkrBaseline=YwkrModel;
GDPBaseline_Young=GDPYoungModel;
EarningsBaseline=EarningsModel;
GDPwkrBaseline=YwkrModel;
LFPBaseline=LFPModel;
WageGapBaseline=WageGapModel;
EarningsBaseline_g=EarningsModel_g;
ConsumpYoungBaseline=ConsumpYoungModel;
EarningsYoungBaseline=EarningsYoungModel;
GDPYoungBaseline=GDPYoungModel;
UtilityBaseline=Utility;
save(['SolveEqmBasic_' CaseName],'GDPBaseline','GDPMktBaseline','GDPwkrBaseline','GDPBaseline_Young','EarningsBaseline','EarningsBaseline_g',...
     'ConsumpYoungBaseline','EarningsYoungBaseline','GDPYoungBaseline','LFPBaseline','WageGapBaseline','UtilityBaseline',...
     'wModel','HModel','HModelAll','pModel');


% % For Benchmark case, let's compute the gains from eliminating TauH, TauW and Both
if isequal(CaseName,'Benchmark');

     [Y_NoTauH,a1,a2,a3,a4,a5,a6,YYoung_NoTauH]=SolveForEqm(zeros(size(TauH)),TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar);    
     [Y_NoTauW,a1,a2,a3,a4,a5,a6,YYoung_NoTauW]=SolveForEqm(TauH,zeros(size(TauW)),Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar);    
     [Y_NoTaus,a1,a2,a3,a4,a5,a6,YYoung_NoTaus]=SolveForEqm(zeros(size(TauH)),zeros(size(TauW)),Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar);    

     Gain_NoTauH=Y_NoTauH./YModel-1;
     Gain_NoTauW=Y_NoTauW./YModel-1;
     Gain_NoTaus=Y_NoTaus./YModel-1;
     YGain_NoTauH=YYoung_NoTauH./GDPYoungModel-1;
     YGain_NoTauW=YYoung_NoTauW./GDPYoungModel-1;
     YGain_NoTaus=YYoung_NoTaus./GDPYoungModel-1;
   
     disp ' '; disp ' ';
     disp '********************************************************************';
     disp '       Counterfactual results with ZERO TAUW AND TAUH';
     disp '********************************************************************'; disp ' ';
     disp ' '; disp ' ';
   
     disp 'OVERALL: Additional output gain over baseline with no frictions (percent):';
     cshow(' ',[Decades 100*[Gain_NoTauH Gain_NoTauW Gain_NoTaus]],'%6.0f %12.1f','Year NoTauH NoTauW NoTauH/W');

     disp ' ';
     disp 'GDPYoung: Additional output gain over baseline with no frictions (percent):';
     cshow(' ',[Decades 100*[YGain_NoTauH YGain_NoTauW YGain_NoTaus]],'%6.0f %12.1f','Year NoTauH NoTauW NoTauH/W');

     % For Baseline case, also show results if Zero TauW/TauH (for Altonji)
     [YModel,YMktModel,EarningsModel,YwkrModel,LFPModel,ConsumpYoungModel,EarningsYoungModel,GDPYoungModel,WageGapModel,EarningsModel_g,Utility,wModel,HModel,HModelAll,pModel,ExitFlag]=SolveForEqm(zeros(size(TauH)),zeros(size(TauW)),Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar);
     SolveEqmBasic_Display

     % Erik's request: display Figure 1 of paper for pModel (zero tauw/tauh)
     % Now we make the figures that Erik wants at the start of the paper

     % Relative propensity
     relpM=zeros(size(pModel)); %pModel=zeros(Noccs,Ngroups,3,Nyears); % YMO is 3rd dimension

     for g=1:Ngroups;
         relpM(:,g,:,:)=pModel(:,g,:,:)./pModel(:,1,:,:);
     end;

     for t=1:Nyears;
     
         % std of log(relp)
         lnrelp_young(:,:,t)=log(relpM(:,:,1,t));
         lnrelp_young(isinf(lnrelp_young))=NaN;
         lnrelp_young(isnan(lnrelp_young))=log(0.001); %min(nanmin(lnrelp_young(:,:,t)));
         meanrelp=nansum(mult(lnrelp_young(:,:,t),earningsweights_avg));
         diff=chadminus(lnrelp_young(:,:,t),meanrelp);
         std_relp_youngZ(t,:)=sqrt(nansum(mult(diff.^2,earningsweights_avg)));  
     end;

     % stdev of relative propensities
     sfigure(3); figsetup;
     greygrid(Decades,(0.5:.5:2));
     plot(Decades,std_relp_young(:,WW),'Color',myblue,'LineWidth',LW); hold on;
     plot(Decades,std_relp_youngZ(:,WW),'Color',mygreen,'LineWidth',LW); hold on;
     %plot(Decades,std_relp_youngZ(:,BM),'Color',mygreen,'LineWidth',LW); 
     %plot(Decades,std_relp_youngZ(:,BW),'Color',mypurp,'LineWidth',LW); 
     vals=(0:.5:2.5)';
     relabelaxis(vals,num2str(vals));
     vals=(1960:10:2010)';
     relabelaxis(vals, num2str(vals),'x');
     chadfig2('  ','  ',1,0);
     if WideFigures; makefigwide; end;
     if 1; %HighQualityFigures;
           %text(1996,1.35,'White women'); %,'Color',myblue);
           %text(1980,0.75,'Black men'); %,'Color',mygreen);
           %text(1967,2.3,'Black women'); %,'Color',mypurp);
         text(1996,.7,'Data');
         text(1980,.3,'Model with zero \tau^w and \tau^h');
     end;
     print('-dpng',['StdlnrelpZeroTaus_' CaseName]);
   
     disp '********************************************************************';
     disp '    END of Counterfactual results with ZERO TAUW AND TAUH';
     disp '********************************************************************';

end;

diary off;
