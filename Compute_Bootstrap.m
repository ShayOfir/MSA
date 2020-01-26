function Bset = Compute_Bootstrap (datum, SV, pdepth, nBS, alpha, varargin)
% Compute confidence interval for Shapley vector using modified bootstrap
% method (with hump analysis)
%
%
%       Bset = Compute_Bootstrap (datum, SV, pdepth, Nboot, alpha, normalize)

if isempty(varargin)
    normalize = 1;
else
    normalize = varargin{1};
end
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
%rng(3); %<--- DEBUG
for qqq=1:nBS
    %debug
    %rng(qqq);
    index1=floor(n*rand(n,1))+1;
    Bset.save_patients(qqq,:) = index1;
    if mod(qqq,50)==0
        disp(['Computing BS-Shapley...',int2str(qqq),'/',int2str(nBS)]);           
    end
    %datum=datamatrix(index1,:);
    SHest = Compute_ShapleyVector_Bound (datamatrix(index1,:), pdepth, normalize);    
    SHest1(:,qqq)=SHest(UU,:)';
end

Bset.Bootstraps = SHest1;
SHAPLEYpartestBOOT=mean(SHest1')';

for jh=1:m %A
    [u, v]=sort(SHest1(jh,:));
     uu=max(u(3),min(u(end-2),u));
    Bset.stdest(jh,:)=[std(u), std(uu)];
end

%rng(3);
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
for j=1:m    %E
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
    mumix2(j)=mu25;pmix1(j)=p15;
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
[calibYY, aver, mode1, factor1] = CalibrateShapleyVector(YY);

Bset.CIcalibmixSHAPL=(100/(mm*m))*[calibYY'-tinv(1-alpha/2,n-1)*Bset.stdestmix(:)*factor,...
                calibYY',...
                calibYY'+tinv(1-alpha/2,n-1)*Bset.stdestmix(:)*factor1,...
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
Bset.factorBOOT = factorBOOT;
calibYYBOOT=aver+factorBOOT*(YY-aver);
Bset.NormCalibSV = (100/(mm*m))*calibYYBOOT;
Bset.aver = aver;
Bset.CIcalibmix=(100/(mm*m))*[calibYYBOOT'-tinv(1-alpha/2,n-1)*Bset.stdestmix(:)*factorBOOT,...
            calibYYBOOT',...
            calibYYBOOT'+tinv(1-alpha/2,n-1)*Bset.stdestmix(:)*factorBOOT,...
            ones(m,1)*mm];
%Bset.STDcalibmix = (100/(mm*m))*tcdf(1,n-1)*Bset.stdestmix(:)*factorBOOT;

Bset.CalibBootstraps = zeros(m,nBS);
for k=1:size(Bset.Bootstraps,2)
    Bset.CalibBootstraps(:,k) = (100/(mm*m))*(aver+factorBOOT*(Bset.Bootstraps(:,k)-aver));
end

end
