%MSA demo

%Initialize:
NR = 26; %No. of brain regions
NP = 50; %No. of patients
TOP = 4.5; %Maximal value of simulated behavior score
pdepth = 4; %perturbation depth
NDAM = 8; %maximum no. of lesioned regions
%Ground-truth function: first six regions are fascilitatory, next 3 are
%inhibitor and the rest are irrelevant
VV=[ones(1,6) -0.5*ones(1,3) zeros(1,NR-9)]; 
alpha = 0.05; %type-I error
nBS = 1000;

WCOND=ones(NP,NR); %init
rng(3); %For demonstration uniform numbers between different runs
% Create virtual pseudo-lesions
for i=2:NP
    Z=1+floor(NR*rand(1,NDAM));
    WCOND(i,Z)=rand(1,NDAM);
end

PERF=WCOND*VV'; %Compute ground-truth performance
xy=[100*(1-WCOND), PERF]; %transform simulated data to the usual data format (numbers are percent and convey damage and not activity)

%Demo #1: Compute Shapley value with Leave-One-Out
[SV, Calib, coal, ~, Lset]=PerformMSA (xy, pdepth, -1, alpha);

%Demo #2: Compute Shapley value with Bootstrap
[SV, Calib, coal, Bset]=PerformMSA (xy, pdepth, nBS, alpha);





figure;
%Bootstrapped, calibrated
subplot(2,2,1);
plotSVandGroundTruth(Bset, 1, 2, 0, alpha, VV);
ylabel ('Calibrated Shapley values');
title ('Bootstrap');

%Bootstrapped, raw
subplot(2,2,2);
plotSVandGroundTruth(Bset, 0, 2, 0, alpha, VV);
ylabel ('Raw Shapley values');
title ('Bootstrap');


%Leave-one-out, calibrated
subplot(2,2,3);
plotSVandGroundTruth(Lset, 1, 2, 0, alpha, VV);
ylabel ('Calibrated Shapley values');
title ('Leave-one-out');


%Leave-one-out, raw
subplot(2,2,4);
plotSVandGroundTruth(Lset, 0, 2, 0, alpha, VV);
ylabel ('Raw Shapley values');
title ('Leave-one-out');


function plotSVandGroundTruth (Bset, calib, fdr_flag, large_only_flag, alpha, VV)

[SV, SVci0, SVci, pval, Z, aver] = StatInference(Bset,calib);
SVsd_pos = SVci - SV;
SVsd_neg = SV - SVci0;    
AsteriskFact = 1.2;
Nreg = length(SV);

hiLim = SVci;
loLim = SVci0;
for j=1:Nreg
    if loLim(j) > -0.5
        loLim(j) = -0.5; 
    end
end

Nreg = length(SV);
LW = 1;
xvec = 1:Nreg;
ast_hi = ones(1,Nreg).*NaN;  
ast_lo = ones(1,Nreg).*NaN;
bigsv = ones(1,Nreg).*NaN;
smallsv = ones(1,Nreg).*NaN;
FDRpval = mafdr(pval,'BHFDR',true);


use_pv2 = ones(1,Nreg);
switch fdr_flag
    case 1 %FDR-corredted p-value
        use_pv = FDRpval;
    case 0 %uncorrected p-value
        use_pv = pval;
    case 2 %show both corrected and un-corrected
        use_pv = FDRpval;
        use_pv2 = pval;
end
ast_hi(use_pv<alpha)=(hiLim(use_pv<alpha))*AsteriskFact;
ast_lo(use_pv<alpha)=(loLim(use_pv<alpha))*AsteriskFact*1.5;
bigsv(Z>0) = 1;
smallsv(Z<0) = 1;
circles = ones(1,Nreg)*NaN;
circles(use_pv2<alpha & bigsv == 1)= (hiLim(use_pv2<alpha & bigsv == 1))*AsteriskFact; %circles are used to denote uncorrected p < 0.05
circles(use_pv<alpha & bigsv ==1) = NaN; %avoid showing dobule marking on the same region
circles = circles .* bigsv;
if ~large_only_flag
   ast_lo = ast_lo.*smallsv;
   ast_hi = ast_hi.*bigsv;
else
    ast_hi = ast_hi.*bigsv;
    ast_lo = ones(1,Nreg)*NaN;
end


h=line([0 Nreg+1],[0 0],'LineStyle',':','Color',[0.5 0.5 0.5],'LineWidth',LW);
aa = h.Parent;
hold on
bar (SV,'FaceColor',[68/256, 114/256, 196/256],'LineWidth',LW);
errorbar (xvec,SV,SVsd_pos,SVsd_neg,'k','LineWidth',LW,'LineStyle','none');               
%a = h.Parent;        
plot (xvec,ast_hi,'LineStyle','none','Marker','*','MarkerSize',12,'MarkerEdgeColor','k'); 
plot (xvec,ast_lo,'LineStyle','none','Marker','*','MarkerSize',12,'MarkerEdgeColor','r');
plot (xvec,circles,'LineStyle','none','Marker','o','MarkerSize',10,'MarkerEdgeColor','k');
yvec = aver.*ones(1,Nreg);
plot ([-1, Nreg+1],ones(1,2).* yvec(1),'k','LineWidth',LW,'LineStyle','--');
%disp (sprintf('Calibrated mean SV=%1.3f',yvec(1)));
aa.LineWidth = LW * 1.25;
aa.XLim= [0 Nreg+1];  
line(xvec,100*VV/sum(VV),'LineWidth',LW,'Color','r');
end

function [SV, SVci0, SVci1, pval, Z, aver] = StatInference (BLset,calib)

    x = BLset{end};
    if isfield(x,'LOO')
        %Leave-one-out
        if calib
            SVci0 = x.CIcalib(:,1);
            SVci1 = x.CIcalib(:,3);
            SV = x.CIcalib(:,2);
        else
            SVci0 = x.CI(:,1);
            SVci1 = x.CI(:,3);
            SV = x.CI(:,2);
        end
        pval = x.pvalest;
        Z = x.Zscoreest;
    else
        %Bootstrap
        if calib
            SVci0 = x.CIcalibmixSHAPL(:,1);
            SVci1 = x.CIcalibmixSHAPL(:,3);
            SV = x.CIcalibmixSHAPL(:,2);
        else
            SVci0 = x.CImix(:,1);
            SVci1 = x.CImix(:,3);
            SV = x.CImix(:,2);
        end   
        pval = x.pvalestmix;
        Z = x.Zscoreestmix;
    end
    
    %The calibrated SVs are normalized by default, but the raw are not.
    %Here we normlize the raw SVs to 100.
    if ~calib
        sumSV = sum(SV);
        SV = 100 * SV ./ sumSV;
        SVci0 = 100 * SVci0 ./ sumSV;
        SVci1 = 100 * SVci1 ./ sumSV;
    end
    aver = sum(SV)/length(SV); 
end
