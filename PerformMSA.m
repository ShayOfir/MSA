
function [SV, Calib, coal, Bset, Lset]=PerformMSA (xy, pdepth, nBS, alpha, TOP, varargin)
%% Performs Multi-perturbation Shapley-value Analysis on a dataset
%
%  [SV, Calib, coal, Bset, Lset] = PerformMSA (xy, pdepth, nBS, alpha, TOP, [optimization, alternative_predictor, normalization]] )
%
% INPUT:
%   xy - matrix in size of (n, m+1); whereas 'n' is number of
%   patients (rows) and 'm' is number of brain regions (columns). The last
%   column correspon to the behavioral score. 
%   Brain regions recieve number between 0 and 100 which correspond to
%   damage percent (float/double etc. are possible).
%   Behavioral score as no limits (againg float/double etc. is possible),
%   as long as higher scores corresponds to better performance (less
%   impairment, and 0 is the minimal score.
%
%   pdepth - bounded perturbation depth. integer between 1 and no. of
%   brain-regions. PAY ATTENTION that high 'pdepth' mat take very long to
%   compute and may exceed memory capcity. As a rule of thumb 15% give reasonable results
%
%   nBS -  statistical inference. integer.
%       = 0 : Compute SV only
%       = -1 : Leave-One-Out method
%       = positive number : number of bootstrap resamples. 1000 is the minimum. 
%
%    alpha - type I error level. Default is 0.05
%
%    optmization: 'gpu' : use GPU if exist, else use parallel CPUs
%    [Default]
%                 'par'  : enforce use of parallel CPUs, even if GPU exists
%                 'none' : do not use any optimization (no parallel CPUs,
%                 no GPUs)
%
%    alternative_predictor - optional. You can override the built-in
%    predictor by providing the Matlab machine-learning model here.
%    Any kind of ML model or even non-model structures, as long as they contain
%    the a method named 'predictFcn' which obey the following sytnax:
%
%    y = alternative_predictor.predictFcn (XX);
%    whereas
%        XX - vector of m brain regions (same as size(xy,2)-1). Each is either 0 (inactive) or 1
%        (active)
%        y - predicted behavioral score, in the same units of the original
%        behavioral scores of xy(:,end)
%
%   normalization - (Default = 1). This option controls the scaling of
%   lesion data. The options are:
%       0 - No scaling (NOT RECOMMENDED! Use it only if you apply your
%       scaling before sending it to PerfromMSA.
%       1 - (Default) Each lesion coloumn is divided by its max 
%       2 - lesion columns are smoothed by dividing them by 10 and rounding
%       before dividing by maximal value. This option is under
%       investigation.

%  OUTPUT:
%   NOTE: Shapley values are computed using potentials method, hence SV of
%   a given depth 'p' is a function of all depths from 1 to 'p'. Hence
%   SV(p,:) means "Shapley value of all coalitions with maximal
%   perturbation depth 'p'"
%
%   SV - raw shaply-vectors. matrix of (pdepth, m). SV(i,:) is the
%   Shapley-vector of pdepth i or less. This is useful for exploratory
%   analysis done to select the appropriate pdepth. Otherwise, use the last
%   row SV(end,:).
%
%   Calib - Calibraite Shapley value. Structure with the following fields:
%      .SV - calibrate Shapley-vectors. matrix of (pdepth, m). 
%      .aver, .mode1, .factor1- calibration factors. vectors at length
%                               'pdepth' (each value corresponds to shapley
%                               vector in SV).
%   The formula for calibration is :  CalibratedSV = aver + (factor1 * rawSV - aver)
%   
%   coal - This is a debug/QC structure and may not be of interest for most
%   users. It allows for following the values generate in the process. Here
%   each cell corresponds only for the associated depth and is not a
%   cumulative measure.
%   The structrure that contains the following fields:
%   .Coal - cell array in length 'pdepth'. Describes the coalitions (in
%           term of indices) which were used in each depth. This is the product of
%           nchosek(1:m,p) for p in 1:pdepth.
%   .Vest - cell array in length 'pdepth'. Describes the predicted performance (V) for each coalition
%           in .Coal
%   .Dist - parameter of the inherent k-NN-based predictor. matrix that decribes the Hamming
%   distance between each coalition in .Coal (no. of rows) to each patient
%   in the dataset (no. of columns). dataset = xy, plus two 'pseudopatients' rows that represent fully-damaged
%   and fully-intact theoretical patietns
%
%  Bset - Bootstrap statistical inference set (produced if nBS>0). This is a cell-array but
%       only the last cell is populated. Always refer to Bset{end}. 
%       Many fields are only relevant for QC / debugging, here are the important
%       fields:
%
%       stdestmix: SD of estimation. Matrix of (m, 2). 1st column is
%       SD. 2nd column is SD without 2 largest and smallest outliers (we
%       recommend using the 2nd column). All of the following computations are
%       done usig the 2nd column!
%
%       zscoreestmix: Z-scores of the shapley-values.
%
%       pvalestmix: p-values of the shapley-values
%
%       pvalestFDRmix: FDR-corrected p-values of the Shapley values.
%
%       CIcalibmixSHAPL: calibrated and normalized (to 100) confidence interval and Shapley-values in an easy to plot version.
%           This is a matrix with size of (m, 4). rows correspons to brain-regions. 
%           1st column: calibrated lower confidence limit(alpha/2)
%           2nd column: calibrated Shapley values (the same as Calib.SV(end,:)
%           3rd column: calibrated higher confidence limit (1-alpha/2)
%           4th column: average (this is the values which SV is calibrated to)
%
%       CImix: the same as CIcalibmixSHAPL but with raw values instead of
%           calibrated
%  
%       Bootstraps, CalibBootstraps - matrix of all SV produced by bootstraping
% 
%       save_patients - matrix of (n, nBS). Counts how many times each patients
%       in the 'xy' was resampled. Values usually span between 0 to 3.
%
%
%  Lset - Leave-one-out statisical inference set (produced in nBS==-1). This is a cell-array but
%         only the last cell is populated. Always refer to Lset{end}. 
%           Many any fields are only relevant for QC / debugging, here are the important
%           fields:
%
%           stdest: SD of estimation. Matrix of (m, 2). 1st column is
%           SD. 2nd column is SD without 2 largest and smallest outliers (we
%           recommend using the 2nd column). All of the following computations are
%           done usig the 2nd column!
%
%           zscoreest: Z-scores of the shapley-values.
%
%           pvalest: p-values of the shapley-values
%
%           pvalestFDR: FDR-corrected p-values of the Shapley values.
%
%           CIcalib: calibrated and normalized (to 100) confidence interval and Shapley-values in an easy to plot version.
%               This is a matrix with size of (m, 3). rows correspons to brain-regions. 
%               1st column: calibrated lower confidence limit(alpha/2)
%               2nd column: calibrated Shapley values (the same as Calib.SV(end,:)
%               3rd column: calibrated higher confidence limit (1-alpha/2)
%       
%
%           CI: the same as CIcali but with raw values instead of
%                   calibrated
%  
%           LOO - matrix of all SV produced by beave-one-outs (the first
%           row correspons to leave out 1st patient and so on)
% 
% COMPATABILITY: The code uses parallel computing toolbox for running the core
% functions on multiple CPU cores. The code was tested on computers with
% intel(R) i5 and i7 with 4, 6 and 8 cores with Matlab versions: R2019b,
% R2020a. The GPU option was tested on NVIDIA Quadro T1000 card with 4GB of memory

%
% Download most updated version from:
% https://github.com/ShayOfir/MSA
%
% Please report any bug to shinofir@gmail.com
%  
% BSD 3-Clause License. Copyright (c) 2020, ShayOfir. All rights reserved.
% For details: https://github.com/ShayOfir/MSA/blob/master/LICENSE.txt

    global nc_data
    nc_data = [];
    
    %GPU Memory parameters
    global chunk_size
    chunk_size = NaN;
    global nchunks
    nchunks = 1;
    normalize_op = 1;
    
    if gpuDeviceCount > 0
            def_optimize = 'gpu';
            [chunk_size, nchunks] = perform_GPU_memory_check(pdepth,xy);
    else
            def_optimize = 'par';
    end

    if isempty(varargin) 
        alt_pred = [];
        optimize = def_optimize;
        normalize_op = 1;
    else
        if length(varargin)>=1
            optimize = varargin{1};
            alt_pred = [];
            if isempty(optimize) || ~(strcmpi(optimize,'gpu') || strcmpi(optimize,'par') || strcmpi(optimize,'none'))
               optimize = def_optimize;
            end
        end
        if length(varargin)>=2
           alt_pred = varargin{2} ;
        end
        if length(varargin)>=3
            normalize_op = varargin{3};
        end
                          
    end
    
    [SV, coal, Calib] = Compute_ShapleyVector_Bound (xy, pdepth(end), TOP, optimize, alt_pred, normalize_op);
    Bset = cell(1,pdepth(end));
    Lset = cell(1,pdepth(end));
    if nBS > 0
        for p=pdepth
            %disp(sprintf('Computing depth %d',p));
            Bset{p} = Compute_Bootstrap(xy, SV, p, nBS, alpha, TOP, optimize, alt_pred, normalize_op);
        end
    end
    if nBS == -1
        for p=pdepth           
            Lset{p} = Compute_LOO(xy, SV, p, alpha, TOP, optimize, alt_pred, normalize_op);
        end        
    end
end

 function [calibYY, aver, mode1, factor1, calib_stat] = CalibrateShapleyVector (YY)
 %% Calibrates Shapley vector
 % Use:
 %
 %  [calibYY, aver, mode1, factor1, calib_stat] = CalibrateShapleyVector (YY)
 %
 % INPUT:
 %        YY - Raw Shapley vector
 %
 % OUTPUT:
 %        calibYY - calibrated Shapley vector
 %        aver - average Shapley value
 %        mode1 - mode.
 %        factor1 - calibration factor
 %        calib_stat: 0= Error, calibration could not be performed
 %                    1= Calibration was successful
 %
 %  The formula for calibration is :  CalibratedSV = aver + (factor1 * rawSV - aver)
  
 mm1g=min(YY);
 MM1g=max(YY);
 perc = 0;
 mode1 = NaN;
 
sigmag=[.125 .25 .50 1]*((MM1g-mm1g)/2);
if sigmag == 0
    %vector is alrady calibtrated
    calibYY = YY;
    aver = mean(YY);
    mode1 = YY(1);
    factor1 = 1;
    calib_stat = 1;
    return
end

xg = zeros(1000,4);
fg = zeros(1000,4);

for j=1:4
    sigg=sigmag(j);
    for i=1:1000
       xg(i,j)=mm1g+(MM1g-mm1g)*i/1000;
       fg(i,j)=sum(exp(-(YY-xg(i,j)).^2/(2*sigg^2)))/sigg;
    end
end
aver=mean(YY);
gg=mean(fg')';
gg4=gg(xg(:,4)<aver);
dgg4=diff(gg4);
[ug4, vg4]=min(dgg4);
if ug4<0 && vg4>0
    [ug, vg]=max(gg4(1:vg4));
    % i=2;while gg(i)>gg(i-1),i=i+1;end,vg=i-1;
    if vg==vg4 || vg==1
        factor1=1;
        calib_stat = 0;
        %disp('no calibration')
    end
    if vg>1 && vg<vg4
        mode1=mm1g+(MM1g-mm1g)*vg/1000;
        factor1=1/(1-mode1/aver);
        calib_stat = 1;
    end
else   
    factor1=1;  
    %disp('no calibration')
    calib_stat = 0;
end
calibYY=(aver+factor1*(YY-aver))*(1+perc*(100/length(YY))/aver);
 end


function Lset = Compute_LOO (datum, SV, pdepth ,alpha, TOP, optimize, alt_pred, normalize_op)
%% Compute confidence interval for Shapley vector using leave-one-out method
%
%       Lset = Compute_LOO (datum, SV, pdepth, alpha [,optmize, alt_pred])
%
% INPUT:
%    datum - same 'xy' in Perform_MSA()
%    SV - raw Shapley value matrix, produced by Compute_ShapleyVector_Bound()
%    pdepth - perturbation depth
%    alpha - type-I error (float/double between 0 and 1)
%    alt_pred - alternative ML predictor. See Perform_MSA() for more
%    information
%
%OUTPUT:
%  Lset - Leave-one-out statistical inference set . This is a cell-array but
%         only the last cell is populated. Always refer to Lset{end}.  
%           Not all fields are relevant for QC / debugging, here are the important
%           fields:
%
%           stdest: SD of estimation. Matrix of (m, 2). 1st column is
%           SD. 2nd column is SD without 2 largest and smallest outliers (we
%           recommend using the 2nd column). All of the following computations are
%           done usig the 2nd column!
%
%           zscoreest: Z-scores of the shapley-values.
%
%           pvalest: p-values of the shapley-values
%
%           pvalestFDR: FDR-corrected p-values of the Shapley values.
%
%           CIcalib: calibrated and normalized (to 100) confidence interval and Shapley-values in an easy to plot version.
%               This is a matrix with size of (m, 3). rows correspons to brain-regions. 
%               1st column: calibrated lower confidence limit(alpha/2)
%               2nd column: calibrated Shapley values (the same as Calib.SV(end,:)
%               3rd column: calibrated higher confidence limit (1-alpha/2)
%       
%
%           CI: the same as CIcali but with raw values instead of
%                   calibrated
%  
%           LOO - matrix of all SV produced by leave-one-outs (the first
%           row correspons to leave out 1st patient and so on)
%


UU = pdepth;
Nbig = 100000;
datamatrix = datum;
[n,m]=size(datamatrix);m=m-1;  
SHAPLEYfullest=SV; 
SHest = zeros(m,n);
Lset.stdest  = zeros(m,2);
II=(1:n)'; 
Lset.save_patients = zeros(n);
SHest1LOO = zeros(m,n);

for qqq=1:n
    disp(['Computing LOO-Shapley...',int2str(qqq),'/',int2str(n)]);           
    datum=datamatrix(II~=qqq,:); %each time leave one patient out
    SHest = Compute_ShapleyVector_Bound (datum, pdepth, TOP, optimize, alt_pred, normalize_op); %compute SV for the leave-one-out dataset
    SHest1LOO(:,qqq) = SHest(UU,:)';    
end
Lset.LOO = SHest1LOO;
SHAPLEYpartestLOO=mean(SHest1LOO')';
%
for jh=1:m
    [u, v]=sort(SHest1LOO(jh,:));
    uu=max(u(3),min(u(end-2),u));
    Lset.stdest(jh,:)=[std(u), std(uu)];
end
Z=randn(n,Nbig);
dfg = zeros(Nbig,2);
for is=1:Nbig
    [u, v]=sort(Z(:,is));
    uu=max(u(3),min(u,u(end-2)));
    dfg(is,:)=[std(u), std(uu)];
end

var=mean(dfg.^2);
Lset.stdest=Lset.stdest*sqrt(n-1)/sqrt(var(2));
mm=mean(SHAPLEYfullest(UU,:));
for jh=1:m
    Lset.Zscoreest(jh)=(SHAPLEYfullest(UU,jh)-mm)/Lset.stdest(jh,2);
    Lset.pvalest(jh)=2*(1-tcdf(abs(Lset.Zscoreest(jh)),n-1));
end
[u, v]=sort(Lset.pvalest);
uu=u*m./(1:m);
uuu(1)=uu(1);
for i=2:m,uuu(i)=max(uuu(i-1),uu(i));end
for i=1:m,Lset.pvalestFDR(v(i))=uuu(i);end
Lset.pval=[Lset.pvalest' Lset.pvalestFDR'];

Lset.SHAPandSTDLOO=[SHAPLEYfullest(UU,:)' Lset.stdest(:,2)]; 
Lset.CI=[SHAPLEYfullest(UU,:)'-tinv(1-alpha/2,n-1)*Lset.stdest(:,2), SHAPLEYfullest(UU,:)', SHAPLEYfullest(UU,:)'+tinv(1-alpha/2,n-1)*Lset.stdest(:,2), ones(m,1)*mm];
YY=SHAPLEYfullest(UU,:);
[calibYY, aver, mode1, factor1] = CalibrateShapleyVector(YY);
Lset.CIcalib=(100/(mm*m))*[calibYY'-tinv(1-alpha/2,n-1)*Lset.stdest(:,2)*factor1, calibYY' calibYY'+tinv(1-alpha/2,n-1)*Lset.stdest(:,2)*factor1];


end


function Bset = Compute_Bootstrap (datum, SV, pdepth, nBS, alpha, TOP, optimize, alt_pred, normalize_op)
%% Compute confidence interval for Shapley vector using modified bootstrap
% method (with hump analysis)
%
%
%       Bset = Compute_Bootstrap (datum, SV, pdepth, nBS, alpha [, alternative_predictor])
%
% INPUT:
%        datum - dataset. see 'xy' in PerformMSA()
%        SV - Raw shapley matrix produced by Compute_ShapleyVector_Bound()
%        pdepth - perturbation depth. See PerformMSA()
%        nBS - no. of boostraps (positive integer)
%        alternative_predictor - see PerformMSA()
%
% OUTPUT:
%       Bset - Bootstrap statistical inference set. This is a cell-array but
%       only the last cell is populated. Always refer to Bset{end}. 
%       Many fields are only relevant for QC / debugging, here are the important
%       fields:
%
%       stdestmix: SD of estimation. Matrix of (m, 2). 1st column is
%       SD. 2nd column is SD without 2 largest and smallest outliers (we
%       recommend using the 2nd column). All of the following computations are
%       done usig the 2nd column!
%
%       zscoreestmix: Z-scores of the shapley-values.
%
%       pvalestmix: p-values of the shapley-values
%
%       pvalestFDRmix: FDR-corrected p-values of the Shapley values.
%
%       CIcalibmixSHAPL: calibrated and normalized (to 100) confidence interval and Shapley-values in an easy to plot version.
%           This is a matrix with size of (m, 4). rows correspons to brain-regions. 
%           1st column: calibrated lower confidence limit(alpha/2)
%           2nd column: calibrated Shapley values (the same as Calib.SV(end,:)
%           3rd column: calibrated higher confidence limit (1-alpha/2)
%           4th column: average (this is the values which SV is calibrated to)
%
%       CImix: the same as CIcalibmixSHAPL but with raw values instead of
%           calibrated
%  
%       Bootstraps, CalibBootstraps - matrix of all SV produced by bootstraping
% 
%       save_patients - matrix of (n, nBS). Counts how many times each patients
%       in the 'xy' was resampled. Values usually span between 0 to 3.
%


Nbig = 100000;
UU = pdepth;
datamatrix = datum;
[n,m]=size(datamatrix);m=m-1;    
SHAPLEYfullest=SV; 
SHest = zeros(m,n);
Bset.stdest  = zeros(m,2);
dfg = zeros(Nbig,2);
II=(1:n)'; 


Bset.save_patients = zeros(nBS,n);
for qqq=1:nBS
    index1=floor(n*rand(n,1))+1;
    Bset.save_patients(qqq,:) = index1;
    if mod(qqq,10)==0 %Prints every 10 bootstraps, thus saving some display space...
        fprintf('Computing BS-Shapley... %d/%d\n',qqq,nBS);           
    end
    SHest  = Compute_ShapleyVector_Bound (datamatrix(index1,:), pdepth, TOP, optimize, alt_pred, normalize_op);    
    SHest1(:,qqq)=SHest(UU,:)';    
end

Bset.Bootstraps = SHest1;
SHAPLEYpartestBOOT=mean(SHest1')';

for jh=1:m 
    [u, v]=sort(SHest1(jh,:));
     uu=max(u(3),min(u(end-2),u));
    Bset.stdest(jh,:)=[std(u), std(uu)];
end

Z=randn(n,Nbig);
dfg = zeros(Nbig,2);
for is=1:Nbig 
    u=sort(Z(:,is));
    uu=max(u(3),min(u,u(end-2)));
    dfg(is,:)=[std(u), std(uu)];
end

var=mean(dfg.^2);
Bset.stdest=Bset.stdest/sqrt(var(2));
mm=mean(SHAPLEYfullest(UU,:));
for jh=1:m  
    Bset.Zscoreest(jh)=(SHAPLEYfullest(UU,jh)-mm)/Bset.stdest(jh,2);
    Bset.pvalest(jh)=2*(1-tcdf(abs(Bset.Zscoreest(jh)),n-1));
end

[u, v]=sort(Bset.pvalest);
uu=u*m./(1:m);
uuu(1)=uu(1); 
for i=2:m
    uuu(i)=max(uuu(i-1),uu(i));
end
for i=1:m
    Bset.pvalestFDR(v(i))=uuu(i);
end
Bset.pval = [Bset.pvalest', Bset.pvalestFDR'];
Bset.SHAPandSTD=[SHAPLEYfullest(UU,:)', Bset.stdest(:,2)]; 
Bset.CI=[SHAPLEYfullest(UU,:)'-tinv(1-alpha/2,n-1)*Bset.stdest(:,2),...
    SHAPLEYfullest(UU,:)',...
    SHAPLEYfullest(UU,:)'+tinv(1-alpha/2,n-1)*Bset.stdest(:,2),...
    ones(m,1)*mm];
YY=SHAPLEYfullest(UU,:);
[calibYY, aver,mode1, factor] = CalibrateShapleyVector(YY);
% In the following line, the (100/(mm*m) is used to normalize the values to
% 0 - 100 range.
Bset.CIcalib=(100/(mm*m))*[calibYY'-tinv(1-alpha/2,n-1)*Bset.stdest(:,2)*factor,...
                calibYY',...
                calibYY'+tinv(1-alpha/2,n-1)*Bset.stdest(:,2)*factor];
            
% Identification of the hump containing the Shapley value
% Initialize some vectors to accelerate things up:
mumix1 = zeros(1,m);
mumix2 = zeros(1,m);
other = zeros(1,m);
pmix1 = zeros(1,m);
pmix2 = zeros(1,m);
sigmix = zeros(1,m);
valleymix = zeros(1,m);
zdistmix = zeros(1,m);

for j=1:m    
    XQ1=sort(SHest1(j,:));
    XQ=XQ1(6:end-5);
    nXQ=length(XQ);
    M2=sum(XQ.^2);
    mu15=mode1;
    mu25=YY(j);
    p15=1/2;
    sig5=std(XQ)/2;
    P1=zeros(size(XQ));
    P2=zeros(size(XQ));
    N1=zeros(size(XQ));
    N2=zeros(size(XQ));
    for EM=1:1000 %F
        p25=1-p15;
        P1=p15*exp(-(XQ-mu15).^2/sig5^2)./(p15*exp(-(XQ-mu15).^2/sig5^2)+p25*exp(-(XQ-mu25).^2/sig5^2));
        P2=p25*exp(-(XQ-mu25).^2/sig5^2)./(p15*exp(-(XQ-mu15).^2/sig5^2)+p25*exp(-(XQ-mu25).^2/sig5^2));
        N1=sum(P1);
        N2=sum(P2);
        S1=sum(XQ.*P1);
        S2=sum(XQ.*P2);
        p15=N1/nXQ;
        p25=N2/nXQ;
        %mu15=S1/max(1,N1);
        mu25=S2/max(1,N2);
        sig5=sqrt((1/nXQ)*(M2+mu15^2*N1+mu25^2*N2-2*(mu15*S1+mu25*S2)));
        %disp([mu1 mu2 mu3 p1 1-p1-p3 p2 sig])
    end    
    mumix1(j)=mu15;
    mumix2(j)=mu25;
    pmix1(j)=p15;
    pmix2(j)=p25;
    sigmix(j)=sig5;
    [mu12, v12]=sort([mumix1(j), mumix2(j)]);
    pt=[p15, p25];
    p12=pt(v12);
    valleymix(j)=(mumix1(j)+mumix2(j))/2-(sigmix(j)^2/(mu12(2)-mu12(1)))*log(pt(2)/pt(1));
    zdistmix(j)=(mu12(2)-mu12(1))/sigmix(j);
    z1=(valleymix(j)-mu12(1))/sigmix(j);
    z2=(-valleymix(j)+mu12(2))/sigmix(j);
    z1=sign(z1)*min(abs(z1),4);
    z2=sign(z2)*min(abs(z2),4);
    fz1=normcdf(z1)-z1*normpdf(z1);fz2=normcdf(z2)-z2*normpdf(z2);
    Bset.stdestmix(j)=Bset.stdest(j,2);
    other(j)=0;
    if zdistmix(j)>1.5 && YY(j)> valleymix(j)
        Bset.stdestmix(j)=std(XQ(XQ>valleymix(j)))/fz2;
        other(j)=mean(XQ(XQ<=valleymix(j)));
    end
    if zdistmix(j)>1.5 && YY(j)< valleymix(j)
        Bset.stdestmix(j)=std(XQ(XQ<valleymix(j)))/fz1;
        other(j)=mean(XQ(XQ>=valleymix(j)));
    end
end

for jh=1:m 
    Bset.Zscoreestmix(jh)=(SHAPLEYfullest(UU,jh)-mm)/Bset.stdestmix(jh);
    Bset.pvalestmix(jh)=2*(1-tcdf(abs(Bset.Zscoreestmix(jh)),n-1));
end
[u, v]=sort(Bset.pvalestmix);
uu=u*m./(1:m);
uuu(1)=uu(1); 
for i=2:m
    uuu(i)=max(uuu(i-1),uu(i));
end
for i=1:m
    Bset.pvalestFDRmix(v(i))=uuu(i);
end


Bset.pvalmix=[Bset.pvalestmix', Bset.pvalestFDRmix'];
Bset.SHAPandSTDmix=[SHAPLEYfullest(UU,:)', Bset.stdestmix(:)];

Bset.CImix=[SHAPLEYfullest(UU,:)'-tinv(1-alpha/2,n-1)*Bset.stdestmix(:),...
    SHAPLEYfullest(UU,:)',...
    SHAPLEYfullest(UU,:)'+tinv(1-alpha/2,n-1)*Bset.stdestmix(:),...
    ones(m,1)*mm];

YY=SHAPLEYfullest(UU,:);
[calibYY, aver, mode1, factor1] = CalibrateShapleyVector(YY);

Bset.CIcalibmixSHAPL=(100/(mm*m))*[calibYY'-tinv(1-alpha/2,n-1)*Bset.stdestmix(:)*factor,...
                calibYY',...
                calibYY'+tinv(1-alpha/2,n-1)*Bset.stdestmix(:)*factor1,...
                ones(m,1)*mm];

Bset.CalibBootstraps = zeros(m,nBS);
for k=1:size(Bset.Bootstraps,2)
    Bset.CalibBootstraps(:,k) = (100/(mm*m))*(aver+factor*(Bset.Bootstraps(:,k)-aver));
end



end

function nc = fast_nchoosek(m,nR,max_nR)
global nc_data
if isempty(nc_data)
    %generate using 'slow' nchoosek:
    for k=1:max_nR
        nc_data{k} = nchoosek(1:m,k);
    end
end
%nchoosek was already calculated, so load the data
nc = nc_data{nR};
end


function [SV, SaveCoal, Calib] = Compute_ShapleyVector_Bound(datum, pdepth, TOP, optimize, alt_pred, normalize_op)
%
%  [SV, SaveCoal, Calib] = Compute_ShapleyVector_Bound(datum, pdepth, vectorize, alternative_predictor, normalize)
%
% Input:
%       datum - matrix of [patients * ROIS; Behavior]
%       pdepth: perturbation depth
%       alternative_predictor: trained ML model object. See PerformMSA().
%

% Output:
%       SV - Shapley Vector
%       SaveCoal - a matrix of working-states
%       Calib - Calibrated Shapley struct:
%           .SV - calibrate Shapley-vectors. matrix of (pdepth, m). 
%           .aver, .mode1, .factor1- calibration factors. vectors at length
%                               'pdepth' (each value corresponds to shapley
%                               vector in SV).
%   The formula for calibration is :  CalibratedSV = aver + (factor1 * rawSV - aver)

if gpuDeviceCount > 0 && strcmpi(optimize,'gpu') && isempty(alt_pred)
    prepare_predictor = 1;
else
    prepare_predictor = 0;
end
if strcmpi(optimize,'none')
    parallel = 0;
else
    parallel = 1;
end

%GPU memory parameters, relevant for GPU only
global chunk_size
global nchunks

datum=Prepare_Dataset_ForPrediction(datum, TOP, normalize_op);
Xdat = datum(:,1:end-1);
[n,m]=size(Xdat);

%The following lines are parameters of the predictor, which were moved here
%for code acceleration  for CPU
BB = 5;
CC = 0;
fast_ones = ones(1,m);
fast_param.CCmat = ones(n,1).*CC;
fast_param.BBmat = ones(n,1).*BB;
fast_param.preXXmat = ones(n,m);

Vdat = datum(:,end);
TOP = max(Vdat);
MAX_DMG_REGIONS = pdepth;
 Vest = cell(1,MAX_DMG_REGIONS);  
if prepare_predictor     
   [Vest, coalitions, SaveCoal] = prepare_predictions(datum, pdepth, chunk_size, nchunks);
else
    %Compute prediction using CPU
    
    % Input number NP of permutations wanted
    
     %alldist = cell(1,MAX_DMG_REGIONS);
     coalitions = cell(1,MAX_DMG_REGIONS);
     for nR = 1:MAX_DMG_REGIONS
         %disp(['Calculating Depth ' num2str(nR)]);
         coalitions{nR} = fast_nchoosek(m,nR,MAX_DMG_REGIONS); %nchoosek(1:m,nR); %All possible coalitions with nRegions out of m
         NCoal = size(coalitions{nR},1);     
         Vest{nR} = zeros(1,NCoal);
         %alldist{nR} = zeros(NCoal,size(Xdat,1));
         %Compute prediction for all coalitions:
         VestPar = Vest{nR};
         coalitionsPar = coalitions{nR};
         DistPar = zeros(size(Xdat,1),NCoal);
         %fprintf('\nApplying predictor for depth=%d',nR);     
         if isempty(alt_pred)
             %K-NN related Predictor
             if parallel
                parfor k=1:NCoal
                    XX = fast_ones;%ones(1,m);
                    XX(coalitionsPar(k,:)) = 0;%Zero the chosen regions  
                    [VestPar(k), DistPar(:,k)] = ApplyPredictor(XX,Xdat,Vdat,fast_param);          
                end
             else
               for k=1:NCoal
                    XX = fast_ones;%ones(1,m);
                    XX(coalitionsPar(k,:)) = 0;%Zero the chosen regions  
                    [VestPar(k), DistPar(:,k)] = ApplyPredictor(XX,Xdat,Vdat,fast_param);          
               end
             end
         else
            for k=1:NCoal
                %Alternative Predictor: trained on damage and not on activity,
                %so XX is reversed
                XX = ones(1,m);
                XX(coalitionsPar(k,:)) = 0;%Zero the chosen regions    
                VestPar(k) = alt_pred.predictFcn(XX);                
            end         
        end
         Vest{nR} = VestPar;
         SaveCoal.Coal{nR} = coalitionsPar;
         SaveCoal.Vest{nR} = VestPar;
         SaveCoal.Dist{nR} = DistPar';
     end   
end

 VVest = zeros(MAX_DMG_REGIONS,m); 
 pre1SHest = zeros(MAX_DMG_REGIONS,m);
 for region=1:m
     for nR=1:MAX_DMG_REGIONS
         %find all coalitions in which region m is present
         if nR < m
            indices = find(sum(coalitions{nR}~=region,2)==nR);
            %Compute mean of(Tij/k-1) 
            VVest(nR,region) = mean(Vest{nR}(indices));
         else
            VVest(nR,region) = 0; %by definition, v(empty coalition) is 0
         end 
            
         
         
         % Compute Sum[Tij/1 + Tij/2 + ... Tij/k]
         if nR == 1 
            pre1SHest(nR,region) = VVest(nR,region)/nR;
         %   preSHest(nR,region) = pre1SHest(nR,region) - VVest(nR,region)/m;
         else  
            pre1SHest(nR,region) = pre1SHest(nR-1,region) + VVest(nR,region)/nR;
         %  preSHest(nR,region) = pre1SHest(nR,region)- VVest(nR,region)/m;
         end                       
     end    
     pre1SHest(MAX_DMG_REGIONS,region) = pre1SHest(MAX_DMG_REGIONS,region)- VVest(MAX_DMG_REGIONS,region)/m;
 end
 %Compute mean TTi
 meanest = zeros(1,MAX_DMG_REGIONS);    
 SHest = zeros(MAX_DMG_REGIONS,m);    
 for j=1:MAX_DMG_REGIONS
    meanest(j)=mean(pre1SHest(j,:));       
    SHest(j,:)=TOP/m+pre1SHest(j,:)-meanest(j);        
 end    
 SV = SHest;
 for j=1:size(SV,1)
    [Calib.SV(j,:), Calib.aver(j), Calib.mode1(j), Calib.factor1(j)] = CalibrateShapleyVector(SV(j,:));
 end
end

function xp = Prepare_Dataset_ForPrediction (x, TOP, varargin)
%% Prepares raw data for prediction: normalizes damage, transform it to extent of activity and adds intact and zero patients
% 
% USE:
%
%      xp = PrepareDatasetForPrediction(x, normalize)
%
% INPUT:
%       x: matrix of n patients X (m regions + 1 behavioral score)
%       normalize:
%                 0= Do no normalize
%                 1= normalize to max(default)
%                 2= smooth and than normalize
%       xp - matrix of n+2 pateints X (m regions + 1 behavioral score) 
% 
% NOTE: We highly recommend normalizing the damage. Our trials showed that
% un-normalized data provides non-informative SV.

[n, m] = size(x);
m = m - 1; 
if isempty(varargin) 
    normalize = 1;
else
    normalize = varargin{1};
end
if normalize > 1
    %smooth before normalizing
    normx = round(x(:,1:end-1) / 10);
    x = [normx, x(:,end)];
end
if normalize > 0
    x=[x(:,1:end-1)./(ones(n,1)*max(x(:,1:end-1))), x(:,end)];%Normalizes damage data
end

x(isnan(x)) = 0;
x(:,1:end-1) = 1 - x(:,1:end-1);
%TOP = max(x(:,end));
xp = x;
xp(n+1,:) = [zeros(1,m), 0]; %Add fully-damaged patient
xp(n+2,:) = [ones(1,m), TOP]; %Add intact patient
end

function [V1est, dista1] = ApplyPredictor (XX, Xdat, Vdat, param)
%% Apply inherent predictor to data
% USE:
% 
%       [Vest, dista1] = ApplyPredictor (XX, Xdat, Vdat)
%
% INPUT:
%       XX - binary coalition at length m (no. of regions)
%       Xdat - working state matrix with +2 imaginary patients (perfectly intact
%               and perfectly damaged): [patients, regions]. Use Prepare_Dataset_ForPrediction() to produce such matrix. 
%       Vat - performance vector associated with Xdat
%
% OUTPUT:
%       Vest - predicted behavioral score
%       dista1 - distance between coalition XX to each of the patients in
%       Xdat


%default value for parameters:


b= 15;
%BB= 5;
%CC =0;  
%TOP = max(Vdat);
n = size(Xdat,1);
%CCmat = ones(n,1)*CC;
%BBmat = repmat(BB,n,1);
%XXmat = repmat(XX,n,1);                
XXmat = XX .* param.preXXmat;
dista1 = min((sum((Xdat - XXmat).^2,2) - param.CCmat),param.BBmat);
W = (exp(-b*dista1))';
gam2 = sum(XX)/sum(W*Xdat);
gam3 = min(1,1/max(gam2*W*Xdat));
V1est=gam3*gam2*sum(W.*Vdat');
end

function [chunk_size, nchunks] = perform_GPU_memory_check (pdepth,xy)
    
    safety_zone = 0.1*10^9; %Residual memory after chunk is stroed
    MemCoef = 13;

    sz = size(xy);
    m = sz(2)-1; %no. of regions, elminate behavior column
    n = sz(1)+2; %no. of patients + 2 pseudo patients
    
    %Compute number of coalitions
    bigNNCoal = 0;
    for k=1:pdepth
        bigNNCoal = bigNNCoal + nchoosek(m,k);
    end
    
    if gpuDeviceCount  == 0 
        % No GPU
        chunk_size = bigNNCoal;
        nchunks = 1;
    else
        % There is GPU
        
        GPU = gpuDevice(1);
        memory_left = GPU.AvailableMemory;

        memory_req = MemCoef * bigNNCoal * m * n;
        if memory_left >= memory_req + safety_zone %required memory is lower than available memory with grace of 200MB        
            chunk_size = bigNNCoal;
            nchunks = 1;
        else
            chunk_size = floor((memory_left - safety_zone)/(MemCoef*m*n));
            nchunks = ceil(bigNNCoal / chunk_size);
        end
    end
disp (['[GPU] No. of memory chunkes=',num2str(nchunks)]);
end

function [Vest, coalitions, SaveCoal] = prepare_predictions (xy, pdepth, chunk_size, nchunks)

%Prepare data for prediction
%pxy = salone_Prepare_Dataset_ForPrediction2(xy);
pxy=xy;
px = pxy(:,1:end-1);
py = pxy(:,end);
[n,m]=size(px);

%Compute all the working-state for the desired perturbation depth
coalitions = cell(1,pdepth);
NCoals = zeros(1,pdepth);

for nR=1:pdepth
    coalitions{nR} = fast_nchoosek(m,nR,pdepth);
    NCoal(nR) = size(coalitions{nR},1);
end

bigNNCoal = sum(NCoal);
bigXX = ones(m,bigNNCoal); %all working states

from_index = 1;
for nR=1:pdepth
    to_index = from_index + NCoal(nR) - 1;
    xx = ones(m,NCoal(nR));
    offset = repmat([0:m:m*NCoal(nR)-1]',1,nR);
    coal_linear = coalitions{nR} + offset;
    xx(coal_linear) = 0;
    bigXX(:,from_index:to_index) = xx;
    from_index = to_index+1;
end

b = 15;
BB=5;
%TOP = max(py);
sz = size(bigXX);
big_aXX = reshape(bigXX,[1,sz(1),sz(2)]); %3D vesion of XX

V1est = zeros(bigNNCoal,1);
from_index = 1;
gpud = struct();
for chunk=1:nchunks
    if chunk < nchunks
        to_index = from_index + chunk_size-1;
    else
        to_index = bigNNCoal;
    end
    NNCoal = to_index - from_index+1;
    XX = bigXX(:,from_index:to_index);
    aXX = big_aXX(:,:,from_index:to_index);
     
    gpud.aXX = gpuArray(single(aXX));
    gpud.XX = gpuArray(single(XX));
    gpud.BBmat = BB .* ones(n,1,NNCoal,'single','gpuArray'); %%
    gpud.px = gpuArray(single(px));
    gpud.Xdat = repmat(gpud.px,[1,1,NNCoal]); %% 
    gpud.py = gpuArray(single(py));
    gpud.XXmat = repmat(gpud.aXX,[n,1,1]); %%

    
    gpud.dista1 = squeeze(arrayfun(@min,sum((arrayfun(@minus,gpud.Xdat,gpud.XXmat)).^2,2),gpud.BBmat));
    gpud.W = pagefun(@transpose,pagefun(@exp,-b.*gpud.dista1));   
    gpud.step5 = pagefun(@transpose,gpud.W*gpud.px);
    gpud.gam2 = sum(gpud.XX,1) ./ sum(gpud.step5,1);
    gpud.gam3 = min(ones(1,NNCoal,'single','gpuArray'),1./max(gpud.gam2.*gpud.step5,[],1));
    gpud.V1est = gpud.gam3.*(gpud.gam2.*sum(gpud.W.*repmat(gpud.py',NNCoal,1),2)');
    V1est(from_index:to_index) = gather(gpud.V1est);
    from_index = to_index+1;
end


from_index = 1;
for nR=1:pdepth
    to_index = from_index + NCoal(nR) - 1;
    Vest{nR} = V1est(from_index:to_index);
    SaveCoal.Coal{nR} = coalitions{nR};
    SaveCoal.Vest{nR} = Vest{nR};
    SaveCoal.Dist{nR} = Vest{nR}.*0; %gathering distance from GPU is inefficient
    from_index = to_index+1;
end




end

