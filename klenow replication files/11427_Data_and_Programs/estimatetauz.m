function [tauw,tauh,z,mg,WageGrowthRelative,EarningsHome]=estimatetauz(alpha,w,phi,s,p,wagebar,earnings,tauhat,Tig,eta,beta,theta,dlta,gam,ConstantAlpha,Tbar,GroupExpAdjustment,WageforTauW,WhichWageGrowth,ConstrainTauH,PurgeWageGrowthSelection,IgnoreBrawnyOccupations,NoFrictions2010)
% estimatetauz.m
%
%  Estimate TauW, TauH, and Z
    
    % -- Dimensions when passed --
    % p(i,g,c,t)
    % wagebar(i,g,c,t)
    % tauhat(i,g,t)
    % Tig(i,g,c,t)
    % Tbar(i,g,c)

    global CaseName OccupationNames ShortNames BrawnyOccupations GroupNames Noccs Ngroups Nyears Ncohorts Decades;
    
    % Setup
    mg=zeros(Ngroups,Nyears);
    tauw=zeros(Noccs,Ngroups,Nyears);
    tauh=zeros(Noccs,Ngroups,Nyears);
    z=zeros(Noccs,Ngroups,Nyears);
    wagegrowth=zeros(Noccs,Nyears);
    wagegrowth_ym=zeros(Noccs,Nyears);
    wagegrowth_mo=zeros(Noccs,Nyears);
    wagegrowthWM=zeros(Noccs,Nyears);
    wagegrowthWM_ym=zeros(Noccs,Nyears);
    wagegrowthWM_mo=zeros(Noccs,Nyears);
    WageGrowthRelative=zeros(Noccs,Ngroups,Nyears);  % Relative wage growth (ww/wm) from Y to M
    ExpAdjustmentFactor_ym=zeros(Noccs,Nyears);
    ExpAdjustmentFactor_mo=zeros(Noccs,Nyears);
    denom=zeros(Noccs,Nyears);
    pfitY=zeros(Noccs,Ngroups,Nyears);
    gg=gam*eta^(eta/(1-eta));
    mu=(1/theta)*(1/(1-eta));
    HOME=1; WM=1; WW=2;
    EarningsHome=squeeze(earnings(HOME,:,:,:)); % g x Cohort x Time
    
    for g=1:Ngroups;
        % First, get taste z and mg(c)
        % mg(c) comes from normalizing z(HOME)=1
        for c=1:Nyears;
            Cohort=7-c;
            pig=p(:,g,Cohort,c);
            mg(g,c)=( wagebar(HOME,g,Cohort,c)*(1-s(HOME,c))^(1/beta) / (gg*Tig(HOME,g,Cohort,c)/Tbar(HOME,g,c)) )^(1/mu);
            mg(g,c)=mg(g,c)*pig(HOME)^(-dlta); % To include taste heterogeneity
            z(:,g,c)=1./(1-s(:,c)) .* ( gg*(pig.^dlta.*mg(g,c)).^mu.*Tig(:,g,Cohort,c)./Tbar(:,g,c) ./ wagebar(:,g,Cohort,c) ).^beta;
            wtilde=Tbar(:,g,c).*w(:,c).*s(:,c).^phi(:,c).*((1-s(:,c)).*z(:,g,c)).^((1-eta)/beta)./tauhat(:,g,c);
            pfitY(:,g,c)=wtilde.^theta / mg(g,c);
        end;

        % Now split tauhat into tauw and tauh
        if ~isnan(ConstantAlpha);
            fprintf('Splitting tauhat with a constant alpha in all years. Alpha = %6.3f\n',alpha);
            tauwall=1-tauhat(:,g,:).^(-alpha); % Constant split in all years
            tauhall=(tauhat(:,g,:).*(1-tauwall)).^(1/eta) - 1;
            TooSmallTauH = (tauhall<ConstrainTauH);  % Perhaps we do not allow tauh to be too close to -1...
            if any(TooSmallTauH~=0);
                %disp 'Stopping to check TooSmallTauH...';
                %keyboard;
                tauhall(TooSmallTauH)=ConstrainTauH;
                tauhatblah=tauhat(:,g,:);
                tauwall(TooSmallTauH)=1-(1+tauhall(TooSmallTauH)).^eta./tauhatblah(TooSmallTauH);
            end;
            tauw(:,g,:)=tauwall;
            tauh(:,g,:)=tauhall;

        else;
            %Use alpha for 1960 and then wage growth for later years
            tauw0=1-tauhat(:,g,1).^(-alpha); % 1960
            tauh0=(tauhat(:,g,1).*(1-tauw0)).^(1/eta) - 1;
            % Check for large negative values of TauH and constrain...
            TooSmallTauH = (tauh0<ConstrainTauH);  % Perhaps we do not allow tauh to be too close to -1...
            if any(TooSmallTauH~=0);
                %disp 'Stopping to check TooSmallTauH...';
                %keyboard;
                tauh0(TooSmallTauH)=ConstrainTauH;
                tauw0(TooSmallTauH)=1-(1+tauh0(TooSmallTauH)).^eta./tauhat(TooSmallTauH,g,1);
            end;
            tauw(:,g,1)=tauw0;
            tauh(:,g,1)=tauh0;

            for c=1:(Nyears-1);
                Cohort=7-c;
                
                % Use WageforTauW = Earnings or Earnings/ReportedHours when estimating TauW (whatever passed)
                wagegrowth_ym(:,c)=WageforTauW(:,g,Cohort,c+1)./WageforTauW(:,g,Cohort,c);         % Y->M
                wagegrowthWM_ym(:,c)=WageforTauW(:,WM,Cohort,c+1)./WageforTauW(:,WM,Cohort,c);     % Y->M
                wagegrowth_mo(:,c)=WageforTauW(:,g,Cohort+1,c+1)./WageforTauW(:,g,Cohort+1,c);     % M->O
                wagegrowthWM_mo(:,c)=WageforTauW(:,WM,Cohort+1,c+1)./WageforTauW(:,WM,Cohort+1,c); % M->O

                % Purge wagegrowth of WW (etc) for differential selection arising from rising LFP
                % as the cohort ages from Y->M or M->O. Assume same rise proportionally in all occs
                % See 2017-12-01-PurgeWageGrowthSelection.pdf
                if PurgeWageGrowthSelection;
                    disp 'PurgeWageGrowth is not fully implemented. Early tests gave weird results';
                    disp 'so exiting now in estimatetauz.m...';
                    keyboard;
                    % Y->M
                    lfpgrowthWM=(1-p(HOME,WM,Cohort,c+1))/(1-p(HOME,WM,Cohort,c));
                    lfpgrowth_g=(1-p(HOME,g, Cohort,c+1))/(1-p(HOME,g, Cohort,c));
                    SelectionEffectYM=(lfpgrowth_g/lfpgrowthWM)^((1-dlta)*mu); % mu=1/theta*1/(1-eta)
                    wagegrowth_ym(:,c)=wagegrowth_ym(:,c)*SelectionEffectYM;
                    % M->O
                    lfpgrowthWM=(1-p(HOME,WM,Cohort+1,c+1))/(1-p(HOME,WM,Cohort+1,c));
                    lfpgrowth_g=(1-p(HOME,g, Cohort+1,c+1))/(1-p(HOME,g, Cohort+1,c));
                    SelectionEffectMO=(lfpgrowth_g/lfpgrowthWM)^((1-dlta)*mu); % mu=1/theta*1/(1-eta)
                    wagegrowth_mo(:,c)=wagegrowth_mo(:,c)*SelectionEffectMO;
                    %if g==2; disp 'purging wage growth...'; keyboard; end;
                end;

                
                % Adjust wagegrowth of WW (etc) for differential experience for purposes of estimating TauW
                % See "GroupExpAdjustment(Noccs,G,Cohort,YMO)" in GetTExperience.m
                % Y->M
                TrueExperience_ym=1+GroupExpAdjustment(:,g,Cohort,2).*(Tig(:,g,Cohort,c+1)-1);
                ExpAdjustmentFactor_ym(:,c)=Tig(:,g,Cohort,c+1)./TrueExperience_ym; % i.e. inflate to give them WM experience, as in model
                wagegrowth_ym(:,c)=wagegrowth_ym(:,c).*ExpAdjustmentFactor_ym(:,c);
                % M->O
                TrueExperience_mo=1+GroupExpAdjustment(:,g,Cohort+1,3).*(Tig(:,g,Cohort+1,c+1)-1);
                ExpAdjustmentFactor_mo(:,c)=Tig(:,g,Cohort+1,c+1)./TrueExperience_mo; % i.e. inflate to give them WM experience, as in model
                wagegrowth_mo(:,c)=wagegrowth_mo(:,c).*ExpAdjustmentFactor_mo(:,c);
                
                %if g==2;
                %disp 'Stopping in estimatetauz. GroupExpAdjustment'; keyboard
                %end;
                
                % Y/M/O cohorts -- use the chosen ones for computing wagegrowth
                if WhichWageGrowth==0; % Use average of both Y->M and M->0
                    wagegrowth=(wagegrowth_ym+wagegrowth_mo)/2;
                    wagegrowthWM=(wagegrowthWM_ym+wagegrowthWM_mo)/2;
                elseif WhichWageGrowth==1; % Use Y->M
                    wagegrowth=wagegrowth_ym;
                    wagegrowthWM=wagegrowthWM_ym;
                elseif WhichWageGrowth==2; % Use M->0
                    wagegrowth=wagegrowth_mo;
                    wagegrowthWM=wagegrowthWM_mo;
                end;
                
                
                OneMinusTauW=(1-tauw(:,g,c)).*wagegrowth(:,c) ./ wagegrowthWM(:,c); %denom(:,c);
                if any(isnan(OneMinusTauW)); % E.g. missing data ==> try just AlphaSplit for the new year
                    missing=isnan(OneMinusTauW);
                    OneMinusTauW(missing)=tauhat(missing,g,c+1).^(-alpha);
                end;
                tauh_prelim=(tauhat(:,g,c+1).*OneMinusTauW).^(1/eta) - 1;

                % Check for large negative values of TauH and constrain...
                TooSmallTauH = (tauh_prelim<ConstrainTauH);  % Perhaps we do not allow tauh to be too close to -1...
                if any(TooSmallTauH~=0);
                    %disp 'Stopping to check TooSmallTauH...';
                    %keyboard;
                    tauh_prelim(TooSmallTauH)=ConstrainTauH;
                    OneMinusTauW(TooSmallTauH)=(1+tauh_prelim(TooSmallTauH)).^eta./tauhat(TooSmallTauH,g,c+1);
                end;
                tauw(:,g,c+1)=1-OneMinusTauW;
                tauh(:,g,c+1)=tauh_prelim; 

                
                % Since we know tauw(HOME)=tauh(HOME)=0, use the wagegrowth equation of model
                % to recover wagegrowth in the HOME sector across the life cycle
                if c<6;
                    wagegrowthHYM = w(1,c+1)./w(1,c) .* Tig(1,g,Cohort,c+1)./Tig(1,g,Cohort,c) .* (s(1,c).^(phi(1,c+1)))./s(1,c).^(phi(1,c));
                    EarningsHome(g,Cohort,c+1)=EarningsHome(g,Cohort,c)*wagegrowthHYM;
                end;
                if c<5;
                    wagegrowthHYO = w(1,c+2)./w(1,c) .* Tig(1,g,Cohort,c+2)./Tig(1,g,Cohort,c) .* (s(1,c).^(phi(1,c+2)))./s(1,c).^(phi(1,c));
                    EarningsHome(g,Cohort,c+2)=EarningsHome(g,Cohort,c)*wagegrowthHYO;
                end;
                
                %disp 'estimatetauz.m home...'; keyboard;
                %if c==1 & g==2; disp 'Firefighters...'; keyboard; end;
                %if g==2; keyboard; end;
            end; % loop over time c
            tauw(HOME,:,:)=0; % Since we don't have wage growth in home sector
            tauh(HOME,:,:)=0;
            tauw(:,WM,:)=0;  % Otherwise wage growth generates a tauw
            tauh(:,WM,:)=0;
        end; % how to split tauhat
        
        % Set frictions to zero for BrawnyOccs or NoFrictions2010, as needed
        if IgnoreBrawnyOccupations; % Only for WW b/c BW are allowed to have frictions
            tauw(BrawnyOccupations,WW,:)=0;
            tauh(BrawnyOccupations,WW,:)=0;            
        end;
        if NoFrictions2010;
            tauw(:,:,6)=0;
            tauh(:,:,6)=0;
        end;
        
        impliedalpha=-log(1-tauw(:,g,:))./log(tauhat(:,g,:)); 
        WageGrowthRelative(:,g,:)=wagegrowth./wagegrowthWM;
        
        % Show values for each group
        disp ' '; disp ' '; 
        disp '--------------------------------------------------------------------------------------------------';
        disp(['                          >>>>> ' GroupNames{g} ' <<<<<'])
        disp '--------------------------------------------------------------------------------------------------';
    
        fmt='%10.3f';
        tle='pY pfitY Z TauHat TauH tauh* TauW0 TauW1 wagegrth wagegrWM ExpFactor alpha';
        for c=1:Nyears;
            Cohort=7-c;
            disp ' ';
            disp '======================================';
            fprintf(' YEAR = %4.0f\n',Decades(c));
            disp '======================================';
            disp '  Note: tauh* := (1+tauh)^eta - 1 ';
            tauhstar=(1+tauh(:,g,c)).^eta-1;
            if c<Nyears; tauwfuture=tauw(:,g,c+1); else tauwfuture=zeros(Noccs,1)*NaN; end;
            cshow(ShortNames,[p(:,g,Cohort,c) pfitY(:,g,c) z(:,g,c) tauhat(:,g,c) tauh(:,g,c) tauhstar tauw(:,g,c) tauwfuture wagegrowth(:,c) wagegrowthWM(:,c) ExpAdjustmentFactor_ym(:,c) impliedalpha(:,c)],fmt,tle);
        end;
        disp ' ';
        disp 'mg(g,t):';
        mg
        
        disp ' '; disp ' ';

        %if g==2; keyboard; end;
        
    end; % loop over groups g

    
  
