% MSA FINAL CODE
function PerformMSA_Analysis (prefix,pdepth,nBS,alpha)
%Performs MSA analysis on a dataset
load ('SubAcuteAndChronic_dataset.mat');
%FMT:
 dataset{2} = SA_rightFMT; nm{2} = 'SArightFMT';
 dataset{1} = SA_leftFMT; nm{1} = 'SAleftFMT';
% dataset{3} = C_leftFMT; nm{3} = 'CleftFMT';
% dataset{4} = C_rightFMT; nm{4} = 'CrightFMT';
% %FMA:
% dataset{5} = SA_rightFMA; nm{5} = 'SArightFMA';
% dataset{6} = SA_leftFMA; nm{6} = 'SAleftFMA';
% dataset{7} = C_leftFMA; nm{7} = 'CleftFMA';
% dataset{8} = C_rightFMA; nm{8} = 'CrightFMA';
% %FMBC:
% dataset{9} = SA_rightFMBC; nm{9} = 'SArightFMBC';
% dataset{10} = SA_leftFMBC; nm{10} = 'SAleftFMBC';
% dataset{11} = C_leftFMBC; nm{11} = 'CleftFMBC';
% dataset{12} = C_rightFMBC; nm{12} = 'CrightFMBC';
% dataset{1} = D_leftFMT; nm{1} = 'DleftFMT';
% dataset{2} = D_rightFMT;nm{2} = 'DrightFMT';
% dataset{1} = same_sa_leftFMT; nm{1} = 'SAleftFMT';
% dataset{2} = same_sa_rightFMT;nm{2} = 'SArightFMT';
% dataset{3} = same_c_leftFMT; nm{3} = 'CleftFMT';
% dataset{4} = same_c_rightFMT; nm{4} = 'CrightFMT';




%pdepth = 4;
%alpha = 0.05;
%nBS = 5;
%xy_prepared = Prepare_Dataset_ForPrediction(rightFMT);
for ds=1:length(dataset)
    disp(nm{ds});
    xy = dataset{ds};
    [SV, coal, d] = Compute_ShapleyVector_Bound (xy, pdepth);
    for p=1:pdepth
        Bset{p} = Compute_Bootstrap(xy, SV, p, nBS, alpha);
    end
    save ([prefix '_' nm{ds} '.mat'],'xy','SV','coal','d','Bset','nBS','alpha','pdepth');
end

end


 function [calibYY, aver, factor, calib_stat] = CalibrateShapleyVector (YY)
 nn=length(YY);mm1g=min(YY);MM1g=max(YY);
sigmag=[.125 .25 .50 1]*((MM1g-mm1g)/2);
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
[ug4 vg4]=min(dgg4);
if ug4<0 & vg4>0
[ug vg]=max(gg4(1:vg4));
% i=2;while gg(i)>gg(i-1),i=i+1;end,vg=i-1;
if vg==vg4 | vg==1
    factor1=1;disp('no calibration')
end
if vg>1 & vg<vg4
    mode1=mm1g+(MM1g-mm1g)*vg/1000;
    factor1=1/(1-mode1/aver);
end
else
    factor1=1;disp('no calibration')
end
%if mode1 >= aver,factor=1;disp('no calibration'),mode1=aver;end
%if vg==length(gg(xg(:,4)<aver)),factor=1;disp('no calibration'),end
%if vg==1,factor=1;disp('no calibration'),end

%figure(1),plot(mm+(MM-mm)*(1:1000)'/1000,[g f]),grid,zoom
%disp([aver mode1 aver*length(YY) factor])
calibYY=aver+factor1*(YY-aver);
 end
%  function [calibYY, aver, factor, calib_stat] = CalibrateShapleyVector (YY)
% % Calibrate Shapley vector:
% %
% %          [calibYY,calib_stat] = CalibrateShapleyVector (YY)
% %
% % input:
% %   YY: Shapley vector
% % output:
% %   calibYY: calibrated Shapley vector
% %   calib_stat: boolean. true= succeed; false= fail.
% %
% 
% nn=length(YY);
% mm1g=min(YY);
% MM1g=max(YY);
% sigmag=[.125 .25 .50 1]*((MM1g-mm1g)/2);
% for j=1:4
%     sigg=sigmag(j);
%     for i=1:1000
%        xg(i,j)=mm1g+(MM1g-mm1g)*i/1000;
%        fg(i,j)=sum(exp(-(YY-xg(i,j)).^2/(2*sigg^2)))/sigg;
%     end
% end
% aver=mean(YY);
% gg=mean(fg')';
% [ug, vg]=max(gg(xg(:,4)<aver));
% % i=2;
% % while gg(i)>gg(i-1),
% %     i=i+1;
% % end
% mode1 = aver;
% if isempty(vg) || (vg == length(gg(xg(:,4)<aver)) || vg == 1)
%     factor = 1;
%     calib_stat = false;
% else
%    mode1=mm1g+(MM1g-mm1g)*vg/1000;
%    factor=1/(1-mode1/aver);
%    calib_stat = true;
% end

% 
% 
% if mode1 >= aver,
%     factor=1;
%     %disp('no calibration'),
%     calib_stat = false;
%     mode1=aver;
% else
%     calib_stat = true;
% end
% calibYY=aver+factor*(YY-aver);
% end

function Bset = Compute_Bootstrap (datum, SV, pdepth, nBS, alpha)
% Compute confidence interval for Shapley vector using modified bootstrap
% method (with hump analysis)
%
%
%       LOOset = Compute_Bootstrap (datum, SV, pdepth, Nboot, alpha)
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
    %debug
    rng(qqq);
    index1=floor(n*rand(n,1))+1;
    LOOset.save_patients(qqq,:) = index1;
    disp(['Computing BS-Shapley...',int2str(qqq),'/',int2str(nBS)]);           
    %datum=datamatrix(index1,:);
    SHest = Compute_ShapleyVector_Bound (datamatrix(index1,:), pdepth);    
    SHest1(:,qqq)=SHest(UU,:)';
end

Bset.Bootstraps = SHest1;
SHAPLEYpartestBOOT=mean(SHest1')';

for jh=1:m %A
    [u, v]=sort(SHest1(jh,:));
     uu=max(u(3),min(u(end-2),u));
    Bset.stdest(jh,:)=[std(u), std(uu)];
end

rng(3);
Z=randn(n,Nbig);
dfg = zeros(Nbig,2);
for is=1:Nbig %B
    u=sort(Z(:,is));
    uu=max(u(3),min(u,u(end-2)));
    dfg(is,:)=[std(u), std(uu)];
end

var=mean(dfg.^2);
Bset.stdest=Bset.stdest/sqrt(var(2));
mm=mean(SHAPLEYfullest(UU,:));
for jh=1:m  %C
    Bset.Zscoreest(jh)=(SHAPLEYfullest(UU,jh)-mm)/Bset.stdest(jh,2);%change
    Bset.pvalest(jh)=2*(1-tcdf(abs(Bset.Zscoreest(jh)),n-1));
end

[u, v]=sort(Bset.pvalest);
uu=u*m./(1:m);
uuu(1)=uu(1); %D
for i=2:m,uuu(i)=max(uuu(i-1),uu(i));end
for i=1:m,Bset.pvalestFDR(v(i))=uuu(i);end
Bset.pval = [Bset.pvalest', Bset.pvalestFDR'];
Bset.SHAPandSTD=[SHAPLEYfullest(UU,:)', Bset.stdest(:,2)]; 
Bset.CI=[SHAPLEYfullest(UU,:)'-tinv(1-alpha/2,m)*Bset.stdest(:,2),...
    SHAPLEYfullest(UU,:)',...
    SHAPLEYfullest(UU,:)'+tinv(1-alpha/2,m)*Bset.stdest(:,2),...
    ones(m,1)*mm];
YY=SHAPLEYfullest(UU,:);
[calibYY, aver, factor] = CalibrateShapleyVector(YY);
Bset.CIcalib=(100/(mm*m))*[calibYY'-tinv(1-alpha/2,n-1)*Bset.stdest(:,2)*factor,...
                calibYY',...
                calibYY'+tinv(1-alpha/2,n-1)*Bset.stdest(:,2)*factor];

% Identification of the hump containing the Shapley value
for j=1:m    %E
    XQ1=sort(SHest1(j,:));
    XQ=XQ1(6:end-5);
    nXQ=length(XQ);
    M2=sum(XQ.^2);
    mu15=mean(XQ)-std(XQ);
    mu25=mean(XQ)+std(XQ);
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
        mu15=S1/max(1,N1);
        mu25=S2/max(1,N2);
        sig5=sqrt((1/nXQ)*(M2+mu15^2*N1+mu25^2*N2-2*(mu15*S1+mu25*S2)));
        %disp([mu1 mu2 mu3 p1 1-p1-p3 p2 sig])
    end    
    mumix1(j)=mu15;
    mumix2(j)=mu25;pmix1(j)=p15;
    pmix2(j)=p25;
    sigmix(j)=sig5;
    [mu12, v12]=sort([mumix1(j), mumix2(j)]);
    pt=[p15 p25];
    p12=pt(v12);
    valleymix(j)=(mumix1(j)+mumix2(j))/2-(sigmix(j)^2/(mu12(2)-mu12(1)))*log(pt(2)/pt(1));
    zdistmix(j)=(mu12(2)-mu12(1))/sigmix(j);
    z1=(valleymix(j)-mu12(1))/sigmix(j);
    z2=(-valleymix(j)+mu12(2))/sigmix(j);
    fz1=normcdf(z1)-z1*normpdf(z1);fz2=normcdf(z2)-z2*normpdf(z2);
    Bset.stdestmix(j)=Bset.stdest(j,2);
    other(j)=0;
    if zdistmix(j)>1.5 & YY(j)> valleymix(j)
        Bset.stdestmix(j)=std(XQ(XQ>valleymix(j)))/fz2;
        other(j)=mean(XQ(XQ<=valleymix(j)));
    end
    if zdistmix(j)>1.5 & YY(j)< valleymix(j)
        Bset.stdestmix(j)=std(XQ(XQ<valleymix(j)))/fz1;
        other(j)=mean(XQ(XQ>=valleymix(j)));
    end
end

for jh=1:m %G
    Bset.Zscoreestmix(jh)=(SHAPLEYfullest(UU,jh)-mm)/Bset.stdestmix(jh);
    Bset.pvalestmix(jh)=2*(1-tcdf(abs(Bset.Zscoreestmix(jh)),n-1));
end
[u, v]=sort(Bset.pvalestmix);
uu=u*m./(1:m);
uuu(1)=uu(1); %H
for i=2:m,uuu(i)=max(uuu(i-1),uu(i));end
for i=1:m,Bset.pvalestFDRmix(v(i))=uuu(i);end


Bset.pvalmix=[Bset.pvalestmix', Bset.pvalestFDRmix'];
Bset.SHAPandSTDmix=[SHAPLEYfullest(UU,:)', Bset.stdestmix(:)]; 
Bset.CImix=[SHAPLEYfullest(UU,:)'-tinv(1-alpha/2,n-1)*Bset.stdestmix(:),...
    SHAPLEYfullest(UU,:)',...
    SHAPLEYfullest(UU,:)'+tinv(1-alpha/2,n-1)*Bset.stdestmix(:),...
    ones(m,1)*mm];
YY=SHAPLEYfullest(UU,:);
[calibYY, aver, factor] = CalibrateShapleyVector(YY);
Bset.CIcalibmixSHAPL=(100/(mm*m))*[calibYY'-tinv(1-alpha/2,n-1)*Bset.stdestmix(:)*factor,...
                calibYY',...
                calibYY'+tinv(1-alpha/2,n-1)*Bset.stdestmix(:)*factor,...
                ones(m,1)*mm];
%I
rt2=SHest1(1,:);
for i=2:m,
    rt2=[rt2, SHest1(i,:)];
end
XQ2=rt2(10:end-9);
mu10=mean(XQ2)-std(XQ2);
mu20=mean(XQ2);
mu30=mean(XQ2)+std(XQ2);
p10=1/3;
p30=1/3;
sig0=std(XQ2)/2;
%J
[mu1, mu2, mu3] = ComputeMixture(XQ, mu10, mu20, mu30, p10, p30, sig0);
mumin=min([mu1 mu2 mu3]);
factorBOOT=1/(1-mumin/aver);
calibYYBOOT=aver+factorBOOT*(YY-aver);
Bset.CIcalibmix=(100/(mm*m))*[calibYYBOOT'-tinv(1-alpha/2,n-1)*Bset.stdestmix(:)*factorBOOT,...
            calibYYBOOT',...
            calibYYBOOT'+tinv(1-alpha/2,n-1)*Bset.stdestmix(:)*factorBOOT,...
            ones(m,1)*mm];
Bset.STDcalibmix = (100/(mm*m))*tcdf(1,n-1)*Bset.stdestmix(:)*factorBOOT;
dbug=0;
end

function [mu1, mu2, mu3] = ComputeMixture (XQ, mu10, mu20, mu30, p10, p30, sig0)
%input data XQ
%input initial guess mu10 mu20 m30 p10 p30 sig0
n=length(XQ);
M2=sum(XQ.^2);
mu1=mu10;
mu2=mu20;
mu3=mu30;
pq1=p10;
pq3=p30;
sig=sig0;
P1=zeros(size(XQ));
P2=zeros(size(XQ));
P3=zeros(size(XQ));
N1=zeros(size(XQ));
N2=zeros(size(XQ));
N3=zeros(size(XQ));
for EM=1:1000
    pq2=1-pq1-pq3;
    P1=pq1*exp(-(XQ-mu1).^2/sig^2)./(pq1*exp(-(XQ-mu1).^2/sig^2)+pq2*exp(-(XQ-mu2).^2/sig^2)+pq3*exp(-(XQ-mu3).^2/sig^2));
    P2=pq2*exp(-(XQ-mu2).^2/sig^2)./(pq1*exp(-(XQ-mu1).^2/sig^2)+pq2*exp(-(XQ-mu2).^2/sig^2)+pq3*exp(-(XQ-mu3).^2/sig^2));
    P3=pq3*exp(-(XQ-mu3).^2/sig^2)./(pq1*exp(-(XQ-mu1).^2/sig^2)+pq2*exp(-(XQ-mu2).^2/sig^2)+pq3*exp(-(XQ-mu3).^2/sig^2));
    N1=sum(P1);
    N2=sum(P2);
    N3=sum(P3);
    S1=sum(XQ.*P1);
    S2=sum(XQ.*P2);
    S3=sum(XQ.*P3);
    pq1=N1/n;
    pq2=N2/n;
    pq3=N3/n;
    mu1=S1/max(1,N1);
    mu2=S2/max(1,N2);
    mu3=S3/max(1,N3);
    sig=sqrt((1/n)*(M2+mu1^2*N1+mu2^2*N2+mu3^2*N3-2*(mu1*S1+mu2*S2+mu3*S3)));
    %disp([mu1 mu2 mu3 p1 1-p1-p3 p2 sig])
end
end

function [SV, SaveCoal, Dist] = Compute_ShapleyVector_Bound(datum, pdepth)
%
%  [SHAPLEYest, SaveCoal, Distance] = Compute_ShapleyVector_Bound(datum, pdepth)
%
% Input:
% datum - matrix of [patients * ROIS; Behavior]
% pdepth: perturbation depth.
%

% Output:
% SV - Shapley Vector
% SaveCoal - a matrix of working-states
% Dist - a matrix of distances from patients



datum=Prepare_Dataset_ForPrediction(datum);
www=size(datum);
nn=www(1); %no. of patients
%XYZ=[(datum(1:nn,1:end-1)./(ones(nn,1)*max(datum(1:nn,1:end-1)))), (datum(1:nn,end))];%Normalizes damage data
%XYZ = datum; %
%b=15; BB=5;CC=0; %initialize parameters

% Input data XYZ of fraction damaged. One more column for performance.
% Zero-performance added for up to 3 working regions. Full-performance
% added for no damage.

% Input number NP of permutations wanted

Xdat = datum(:,1:end-1);
[~,m]=size(Xdat);
Vdat = datum(:,end);
TOP = max(Vdat);
%Coalitions method
if pdepth > 7
    MAX_DMG_REGIONS = 7;
else
    MAX_DMG_REGIONS = pdepth;
end     
 Vest = cell(1,MAX_DMG_REGIONS);  
 alldist = cell(1,MAX_DMG_REGIONS);
 coalitions = cell(1,MAX_DMG_REGIONS);
 for nR = 1:MAX_DMG_REGIONS;
     coalitions{nR} = nchoosek(1:m,nR); %All possible coalitions with nRegions out of m
     NCoal = size(coalitions{nR},1);     
     Vest{nR} = zeros(1,NCoal);
     alldist{nR} = zeros(NCoal,size(Xdat,1));
     %Compute prediction for all coalitions:
     for k=1:NCoal
         XX = ones(1,m);
         XX(coalitions{nR}(k,:)) = 0;%Zero the chosen regions
         [Vest{nR}(k),alldist{nR}(k,:)] = ApplyPredictor(XX,Xdat,Vdat);             
     end
 end     
 VVest = zeros(MAX_DMG_REGIONS,m);     
 for region=1:m
     for nR=1:MAX_DMG_REGIONS
         %find all coalitions in which region m is present
         indices = find(sum(coalitions{nR}~=region,2)==nR);
         %Compute mean of(Tij/k-1)          
         VVest(nR,region) = mean(Vest{nR}(indices));
         % Compute Sum[Tij/1 + Tij/2 + ... Tij/k]
         if nR == 1            
            preSHest(nR,region) = VVest(nR,region);
         else            
            preSHest(nR,region) = preSHest(nR-1,region) + VVest(nR,region)/nR;
         end
     end
 end
 %Compute mean TTi
 meanest = zeros(1,MAX_DMG_REGIONS);    
 SHest = zeros(MAX_DMG_REGIONS,m);    
 for j=1:MAX_DMG_REGIONS
    meanest(j)=mean(preSHest(j,:));       
    SHest(j,:)=TOP/m+preSHest(j,:)-meanest(j);        
 end    
 SV = SHest;   
 SaveCoal = coalitions;
 Dist = alldist;   
end %elementary_MSA

function xp = Prepare_Dataset_ForPrediction (x)
% Prepares raw data for prediction: normalizes damage, transform it to extent of activity and adds intact and zero patients
% 
%      xp=PrepareDatasetForPrediction(x)
%
% x - matrix of n patients X (m regions + 1 behavioral score)
% xp - matrix of n+2 pateints X (m regions + 1 behavioral score) 
[n, m] = size(x);
m = m - 1; 
x=[x(:,1:end-1)./(ones(n,1)*max(x(:,1:end-1))), x(1:n,end)];%Normalizes damage data
x(:,1:end-1) = 1 - x(:,1:end-1);
TOP = max(x(:,end));
xp = x;
xp(n+1,:) = [zeros(1,m), 0];
xp(n+2,:) = [ones(1,m), TOP];
end

function [Vest, dista1] = ApplyPredictor (XX, Xdat, Vdat)
%
% [Vest, distance] = ApplyPredictor (XX, Xdat, Vdat)
%
% XX - binary coalition at length m (no. of regions)
% Xdat - working state matrix with +2 imaginary patients (perfectly intact
% and perfectly damaged): [patients, regions]
% Vat - performance vector associated with Xdat
%

%default value for parameters:
%[n,m]=size(Xdat);
b= 15;
BB= 5;
CC =0;
TOP = max(Vdat);
n = size(Xdat,1);
CCmat = ones(n,1)*CC;
BBmat = repmat(BB,n,1);
XXmat = repmat(XX,n,1);                
dista1 = min((sum((Xdat - XXmat).^2,2) - CCmat),BBmat);
W = (exp(-b*dista1))';
WXdat = W * Xdat;
%gam1 = (WXdat * XX') / (WXdat * WXdat');
gam2 = sum(XX)/sum(WXdat);
gam3 = min(1,1/max(gam2*W*Xdat));
sum_W = sum(W); 
bet=1/sum_W;
sumWVdat_t = sum(W.*Vdat'); 
Vest=min(TOP,gam2*sumWVdat_t);
end

