%MSA demo

%Load clinical data
dat = load('ClinicalCohorts.mat');
%The MAT file contains the following variables:
% LHD_xy = matrix of lesion and Fugl-Meyer Assessment score of the upper
% limb of left-hemispheric damaged (LHD) patients 
% RHD_xy = the same for right-hemispheric damaged (RHD) patients
% RegionNames - cell-vector of strings which includes the region names of
% the 26 preselected region-of-interest
%
% The input matrices should be contructed as follows:
% Each row represents a patient. Each column (except of the last) represent
% as brain region (in the same oreder as 'RegionNames'. Number represents
% percents (0 to 100). I used double percision, but other percisions might
% also work.
% The last column is the behavioral score. The convention is that higher
% numbers represent better performance, i.e. less impairment. Use
% transformation to represnts the opposite (e.g. max-score - x)

pdepth = 5;
%Set perturbation depth (J) to 5. This means that only coalitions with 5
%inactive regions or less would be taken into account.
% WARNING: HIGH pdepths may require considerable computation time and may
% exhoust RAM! To get some intuiation use can use nchoosek(K,J) where K is
% the number of regions and J is perturbation depth. For the current
% cohorts K = 26, hence J = 5 yields 83,681 different coalitions while J =
% 10 will yield > 10M coalitions, putting high burden on system resources.
 
nBS = 1000;
%Set no. of bootstrap resamples to 1000. This is the minimal number that
%yield good results. Higher number require linearly longer computation
%time.

alpha = 0.05;
%Set significance level to 0.05

TOP = 66;
%TOP must be the expected score of an intact individual. The FMA spans between 0 and 66.
%Pay attention that %zero is always the expected score a maximally injured patient. In case you
%score minimum is not 0, it is advisiable to use a transformation (i.e.,
%the result will represent the score that patient achieved above the
%minimum)

optimize = 'gpu';
%The code is planned to use GPU whenever possible, and if no GPU is detected it
%uses parralel CPU cores. Only if these two are not possible classical
%non-paralelized code is invoked. Use 'par' to enforce use of parallel CPU
%cores and 'none' to enforce no-optimization

%For the sake of demonstration, We use here only the right-hemispheric cohort. 

%Demo #1: Compute Shapley value with Leave-One-Out (Jacknife method)
[SV, Calib, ~, ~, Lset]=PerformMSA (RHD_xy, pdepth, -1, alpha, TOP, optimize);


%Demo #2: Compute Shapley value with Bootstrap
[SV, Calib, ~, Bset]=PerformMSA (RHD_xy, pdepth, nBS, alpha, TOP, optimize);

%There is some redundency, as Bset was added in subsequent stage. To get
%the Shapley-vector information use can use SV (for raw SV) or Calib.SV
%(for calibrated SV). Bset and Lset include in addition to the raw and calibrated
% SVs also statistical data derived from jacknife/bootstrap methods. See
% the help for PerformMSA.m for details.

save ('msa_demo_data.mat');



figure;
%Bootstrapped, calibrated
subplot(2,2,1);
plotSV(Bset, 1, 2, 1, alpha, dat.RegionNames);
ylabel ('Calibrated Shapley values');
title ('Bootstrap');

%Bootstrapped, raw
subplot(2,2,2);
plotSV(Bset, 0, 2, 1, alpha, dat.RegionNames);
ylabel ('Raw Shapley values');
title ('Bootstrap');


%Leave-one-out, calibrated
subplot(2,2,3);
plotSV(Lset, 1, 2, 1, alpha, dat.RegionNames);
ylabel ('Calibrated Shapley values');
title ('Leave-one-out');


%Leave-one-out, raw
subplot(2,2,4);
plotSV(Lset, 0, 2, 1, alpha, dat.RegionNames);
ylabel ('Raw Shapley values');
title ('Leave-one-out');


function plotSV (Bset, calib, fdr_flag, large_only_flag, alpha, RegionNames)

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
%bar (SV,'FaceColor',[68/256, 114/256, 196/256],'LineWidth',LW);
plot (SV,'Color',[68/256, 114/256, 196/256],'LineWidth',LW);
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
aa.XTick = 1:26;
aa.XTickLabel = RegionNames;
aa.XTickLabelRotation = 90;
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
