% ShowParameters.m
%
%  Show the key parameters for the Talent project

mu=1/theta*1/(1-eta);
GammaBase=gamma(1-1/(theta*(1-eta)));
GammaTilde=exp(0.57721566/(theta*(1-eta)));
GammaBar=GammaTilde^(1-dlta)*GammaBase^dlta;
gam=GammaBar;
   
if FiftyFiftyTauHat;
    AlphaSplitTauW1960=1/2;
end;
disp ' ';
disp '=============================================================';
disp(['KEY PARAMETER VALUES:    CaseName = ' CaseName]);;
fprintf('                   theta =%8.4f\n',theta);
fprintf('                     eta =%8.4f\n',eta);
fprintf('            theta*(1-eta)=%8.4f\n',theta*(1-eta));
fprintf('                    dlta =%8.4f\n',dlta);
fprintf('                      mu =%8.4f\n',mu);
fprintf('                   sigma =%8.4f\n',sigma);
fprintf('           EstimateDelta =%8.4f\n',EstimateDelta);
fprintf('           ConstrainTauH =%8.4f\n',ConstrainTauH);
fprintf('PurgeWageGrowthSelection =  %4.0f\n',PurgeWageGrowthSelection);

if exist('AlphaSplitTauW1960');
    fprintf('      AlphaSplitTauW1960 =%8.4f\n',AlphaSplitTauW1960);
end;
fprintf([  '   OccupationforWageHome = %2.0f ' OccupationNames(OccupationtoIdentifyWageHome,:) '\n'],OccupationtoIdentifyWageHome);
% if exist('PhiKeyOcc');
%     fprintf('       phi(KeyOcc) = '); fprintf('%7.4f',PhiKeyOcc); disp ' ';
% end;
% fprintf('         phi(Home) = '); fprintf('%7.4f',phiHome); disp ' ';
if FiftyFiftyTauHat;
    fprintf('FiftyFiftyTauHat = 1 -- robust 50/50 split of tauhat\n');
end;
if IgnoreBrawnyOccupations;
    fprintf('IgnoreBrawnyOccupations = 1 -- choosing T(i,g) to zero out tauw/h there\n');
end;
if SameExperience==0;
    fprintf('SameExperience = 0 -- all occs have different returns to experience\n');
end;
if ~isnan(ConstantAlpha);
    fprintf('   ConstantAlpha =%8.4f\n',ConstantAlpha);
end;
if ~isnan(Alpha0FixedSplit);
    fprintf('   Alpha0FixedSplit =%8.4f\n',Alpha0FixedSplit);
end;
if NoFrictions2010;
    fprintf('NoFrictions2010 = 1 -- choosing T(i,g,2010) s.t. set tauw(2010)=tauh(2010)=0\n');
end;
if WageGapAdjustmentFactor~=1;
    fprintf('WageGapAdjustmentFactor =%8.4f\n',WageGapAdjustmentFactor);
    if WageGapAdjustmentFactor==0;
        disp('  (i.e. using WM wages for all groups, so no wage gaps');
    end;
end;
if HalfExperience;
    NoGroupAdj=1; % Do not do the other group adjustment
    fprintf('HalfExperience= 1 -- Give WW/BM/BW half the returns to experience');
end;
if NoGroupExpAdj;
    fprintf('NoGroupExpAdj = 1 -- Do not adjust wage growth for experience Y->M\n');
end;
if WhichWageGrowth==1;
    fprintf('WhichWageGrowth==1 -- Use Y->M wage growth when estimating TauW\n');
end;
if WhichWageGrowth==2;
    fprintf('WhichWageGrowth==2 -- Use M->O wage growth when estimating TauW\n');
end;
if isequal(CaseName,'TauWWisZero');
    disp ' '; 
    disp '********************************************************';
    disp 'NOTE WELL: WW and WM are *swtiched* in everything that follows';
    disp '   in order to check robustness to assuming tau(WW)=0 as our';
    disp '   normalization';
    disp '********************************************************';
    disp ' '; 
end;

disp '=============================================================';
disp ' ';
