function [w,H,wtilde,HModelAll,pmodel]=solveeqm(x,t,TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,dlta,mu,sigma,Tbar); 

% function [w,H,wtilde,HModelAll,pmodel]=solveeqm(x,t,TauH,TauW,Z,TExperience,A,phi,q,wH_T,gam,GammaBase,beta,eta,theta,mu,sigma,Tbar); 
%
% Given a guess for x=[mgtilde Y] 5x1 and a year t (e.g. 1=1960), 
% solve for w(i) in year t
%
% Returns w,H ==> Noccs x 1   and wtilde ==> Noccs x Ngroups
%    and HModelAll, pmodel ==> Noccs x Ngroups x YMO
%
% Updated for Version 6 in December 2018 to include dlta / taste heterogeneity

global Noccs Ngroups Ncohorts Nyears CohortConcordance TauW_Orig pData HAllData 
global TauW_C phi_C mgtilde_C w_C GeoArith_adjust StopHere % For keeping track of history in solution

mgtilde_t=x(1:4);
Y_t=x(5);

% Update the History variables with our candidate answers
ct=7-t; % ct := The "cohort" entry (1..8) corresponding to Year t
mgtilde_C(:,ct)=mgtilde_t;
w=zeros(Noccs,1);
H=zeros(Noccs,1);
wtilde=zeros(Noccs,Ngroups);
HModelAll=zeros(Noccs,Ngroups,3); %Noccs x Ngroups x YMO
pmodel=zeros(Noccs,Ngroups,3); %Noccs x Ngroups x YMO
Pwork=zeros(Noccs,Ngroups,3); %Noccs x Ngroups x YMO
                              
%w0=[5000 10000]; fzerofactor=[2 1.4]; NumTries=8;
w0=[1000 4000]; fzerofactor=[2 2]; NumTries=8;
if theta>3;
    w0=[2000 4000]; fzerofactor=[2 1.4]; NumTries=8;
end;
if dlta>0;
    w0=[5000 10000];
end;

%i=1
%e_HSupplyDemand(797)
%e_HSupplyDemand(w0(1))
%e_HSupplyDemand(w0(2))
%keyboard

%wi=fzerochad(@e_HSupplyDemand,w0,fzerofactor,NumTries,1)
% for i=2:Noccs;
%     i
%     e_HSupplyDemand(3850)
% end;
%keyboard
%abc
 
for i=1:Noccs;   %2:Noccs; % Excluding occ=1 Home
    %fprintf('Solving Occupation i = %2.0f\n',i);
    fSupplyDemand=@(w) e_HSupplyDemand(w);
    %if i==25; disp 'i=25 stopping'; keyboard; end;
    wi=fzerochad(fSupplyDemand,w0,fzerofactor,NumTries);
    [resid,Hi,wtildei,HAlli,pmodeli]=e_HSupplyDemand(wi);
    w(i)=wi;
    H(i)=Hi;
    wtilde(i,:)=wtildei;
    pmodel(i,:,:)=pmodeli; 
    HModelAll(i,:,:)=HAlli;
    %if i==23 && StopHere==1;
    %    StopHere=2;
    %    fSupplyDemand(wi);
    %end;
    w0=[.7*wi 2*wi];
end;
%pmodel(1,:,:)=1-sum(pmodel(2:Noccs,:,:)); % fill in Home for LFP  10/25/17--just like any other occ
%if any(any(sum(pmodel(1:Noccs,:,:))~=1));
%    disp 'sum(pmodel(1:N)) ~=1????. Maybe okay since mg(c) not solved out???...'; %keyboard;
%end;

% ----------------------------------
% NESTED FUNCTIONS (inherit variables)
% ----------------------------------

    function [resid,Hdemand,wtilde_i,HAll_i,pmodel_i]=e_HSupplyDemand(wi);  % Nested function -- can see the other variables in play.
    
        Hdemand=A(i,t).^(sigma-1)./wi^sigma*Y_t;
        Hsupply=0;
        w_C(i,ct)=wi; % Update with candidate
                
        % Hsupply requires more work. Also, special treatment for 1960 and 1970 because of 1940/50 cohorts
        % We take the basic data as given for 1950/1940 cohorts and only adjust LF participation bc of TauW
        
        for ymo=0:2; % Loop over Y/M/O cohorts in year t. All groups at same time as 1x4 vectors
            c=CohortConcordance(t,2+ymo); % Cohort index for YMO in year t
            % Pull out the relevant parameters -- date t
            phi_t=phi(i,t);
            tauw_t=TauW(i,:,t);
            
            if c==7 | c==8; 
                % 1940/50 cohort ==> use pig from observed year (e.g. 1960) since not changing 
                % over life cycle in model
                pig_t=squeeze(pData(i,:,c,t));
                HAll_i(:,1+ymo)=HAllData(i,:,c,t); % Data

 % For testZ.m               
 %disp 'Chad testing fix in solveeqm.m. Drop this later!!!!'; 
 %HAll_i(:,1+ymo)=HAll_i(:,ymo);  % Just use Young instead of data
                
                
            else; % No more special cases: c={6,5,4,3,2,1}
                
                % When young ==> cohort c
                phi_c=phi_C(i,c);  % NxC cohort version of phi 
                s_c=1/(1+(1-eta)/beta/phi_c);
                tauh_c=TauH(i,:,c);
                tauw_c=TauW_C(i,:,c); 
                z_c=Z(i,:,c);
                
                tau_c=(1+tauh_c).^eta./(1-tauw_c);
                w_c=w_C(i,c);   % w_C should be NxC = Nx8
                mgtilde_c=mgtilde_C(:,c)';

                % Exogenous experience 
                texp_c=squeeze(TExperience(i,:,c,t-ymo));
                texp_t=squeeze(TExperience(i,:,c,t));
                tbar_c=Tbar(i,:,t-ymo);

                wtilde_c=tbar_c*w_c*s_c^phi_c.*((1-s_c).*z_c).^((1-eta)/beta)./tau_c;                
                mg_c=mgtilde_c.^(1/mu); %IMPORTANT: mgtilde := mg^(1/theta*1/(1-eta)), so mg=mgtilde^(1/mu)
                pig_t=wtilde_c.^theta ./ mg_c;
                
                % Arithmetic mean E[h*e|Work] = term1(c)*term2(t)
                % Updated to include dlta
                term1=GammaBase*(eta*s_c^phi_c*tbar_c*w_c.*(1-tauw_c)./(1+tauh_c)).^(eta/(1-eta));
                %term1=gam*(eta*s_c^phi_c*tbar_c*w_c.*(1-tauw_c)./(1+tauh_c)).^(eta/(1-eta));
                term2=s_c^phi_t*( (1-dlta)*(1./pig_t).^mu + dlta);  % mu:=1/theta*1/(1-eta)
                AvgQuality=term1.*term2;
     
                % Chad Adjustment testing 1/9/19 to adjust for geo/arith in data vs model
                % AvgQuality=.6893/.9*AvgQuality;
                % This is just effectively a *level* normalization in 1960
                % (e.g. we could do it to the *data* for 1960/1970 Middle/Old instead)
                %
                %  Eqn (B5) in the appendix shows that in our model, the ratio of the 
                %  geometric mean to the arithmetic mean is
                %   G = A * gam/GammaBase * pig^(dlta*mu) / ((1-dlta)+dlta*pig^mu);
                %  1/10/19 old line in solveeqm.m:    AvgQuality=AvgQuality*gam/GammaBase/GeoArith_markup;
                %  1/16/19 new line:                  AvgQuality=AvgQuality*GeoArith_adjust;
                   
                AvgQuality=AvgQuality*GeoArith_adjust;
                
     % % Rough Fix 1/10/19 NEED TO FORMALIZE!!!!!!!
     % % THIS IS JUST A SHORTCUT TO SEE HOW THE RESULTS LOOK
     % %  For dlta=1/2, based on EarningsData/EarningsModel for 1990-2010=1.5846
     % if dlta==1/2;
     %     AvgQuality=AvgQuality/1.5846; 
     % end;
                
                HAll_i(:,ymo+1)=(q(:,c,t)'.*pig_t.*texp_t.*AvgQuality)';
                if any(isnan(HAll_i(:,ymo+1))); disp 'isnan(HAll_i)'; keyboard; end;

            end;
                
            % Hsupply
            pmodel_i(:,ymo+1)=pig_t'; % G x YMO
            Hsupply=Hsupply+sum(HAll_i(:,ymo+1));
            if ymo==0; 
                wtilde_i=wtilde_c;
                %Hyoung_i=Hsupply;
            end;
            
            if ~isreal(Hsupply*Hdemand) | isnan(Hsupply*Hdemand);
                disp 'imaginary/NaN Hsupply or Hdemand. stopping in solveeqm function'
                keyboard
            end;
            %if StopHere==2;
            %  disp 'Stopping bc of StopHere...'; keyboard;
            %end;
            
        end; % sum over cohorts c
        resid=Hsupply-Hdemand;
        
    end; % Function e_HSupplyDemand
     
    % Main function w=solveeqm(x,t)  -- nested functions require that main function have "end" 
end