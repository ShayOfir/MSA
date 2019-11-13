% This program accept one damage-performance matrix datamatrix and
% calculates from subjects BOOTSTRAP standard errors and p-values
alpha=.1;%UU=4;
NN1=1000;
clear datum SHAPLEYfullestBOOT SHAPLEYpartestBOOT other II
clear SHest1BOOT u v uu uuu stdestBOOT pvalestBOOT pvalestFDRBOOT  
clear CIBOOT ZscoreestBOOT asd pvalBOOT SHAPandSTDBOOT CIcalibBOOT YY pvalestFDRnewBOOT RT
clear CInewBOOT ZscoreestnewBOOT pvalnewBOOT SHAPandSTDnewBOOT CIcalibnewBOOT stdestnewBOOT pvalestnewBOOT
clear CImixBOOT ZscoreestmixBOOT pvalmixBOOT SHAPandSTDmixBOOT CIcalibmixBOOT stdestmixBOOT pvalestmixBOOT
clear mumix1 mumix2 pix1 pmix2 sigmix valleymix zdistmix pvalestFDRmixBOOT
III=3;
datum=datamatrix;
asd=size(datamatrix);
MSAcoalitions_est
for j=1:UU
SHAPLEYfullestBOOT(j,:)=SHest(j,:);
end
clear datum II,II=(1:asd(1))'; 

    for qqq=1:NN1
%    datum=datamatrix(II~=qqq,:);
index1=floor(asd(1)*rand(asd(1),1))+1;
datum=datamatrix(index1,:);
if 100*floor(qqq/100)==qqq, disp(qqq),end
    MSAcoalitions_est
  for j=1:UU
    SHest1BOOT(:,qqq)=SHest(j,:)';
 end

    SHAPLEYpartestBOOT=mean(SHest1BOOT')';
end
for jh=1:asd(2)-1,
    [u v]=sort(SHest1BOOT(jh,:));
    uu=max(u(3),min(u(end-2),u));
    stdestBOOT(jh,:)=[std(u) std(uu)];
end
    
clear Z u v uu
Z=randn(asd(1),100000);
for is=1:100000
    [u v]=sort(Z(:,is));
    uu=max(u(3),min(u,u(end-2)));
    dfg(is,:)=[std(u) std(uu)];
end
var=mean(dfg.^2);
stdestBOOT=stdestBOOT/sqrt(var(2));

mm=mean(SHAPLEYfullestBOOT(UU,:));
for jh=1:asd(2)-1
    ZscoreestBOOT(jh)=(SHAPLEYfullestBOOT(UU,jh)-mm)/stdestBOOT(jh,2);
    pvalestBOOT(jh)=2*(1-tcdf(abs(ZscoreestBOOT(jh)),asd(1)-1));
end
clear u v uu uuu
[u v]=sort(pvalestBOOT);
uu=u*(asd(2)-1)./(1:asd(2)-1);
uuu(1)=uu(1);
for i=2:asd(2)-1,uuu(i)=max(uuu(i-1),uu(i));end
for i=1:asd(2)-1,pvalestFDRBOOT(v(i))=uuu(i);end

pvalBOOT=[pvalestBOOT' pvalestFDRBOOT'];
SHAPandSTDBOOT=[SHAPLEYfullestBOOT(UU,:)' stdestBOOT(:,2)]; 
CIBOOT=[SHAPLEYfullestBOOT(UU,:)'-tinv(1-alpha/2,asd(1)-1)*stdestBOOT(:,2) SHAPLEYfullestBOOT(UU,:)' SHAPLEYfullestBOOT(UU,:)'+tinv(1-alpha/2,asd(1)-1)*stdestBOOT(:,2) ones(asd(2)-1,1)*mm];

YY=SHAPLEYfullestBOOT(UU,:);gausskervalley
CIcalibBOOT=(100/(mm*(asd(2)-1)))*[calibYY'-tinv(1-alpha/2,asd(1)-1)*stdestBOOT(:,2)*factor1 calibYY' calibYY'+tinv(1-alpha/2,asd(1)-1)*stdestBOOT(:,2)*factor1];

% Identification of the hump containing the Shapley value
for j=1:asd(2)-1
    clear XQ XQ1 P1 P2 N1 N2 S1 S2
XQ1=sort(SHest1BOOT(j,:));XQ=XQ1(6:end-5);
n=length(XQ);M2=sum(XQ.^2);
%mu15=mean(XQ)-std(XQ);mu25=mean(XQ)+std(XQ);p15=1/2;sig5=std(XQ)/2;
mu15=mode1;mu25=YY(j);p15=1/2;sig5=std(XQ)/2;
P1=zeros(size(XQ));P2=zeros(size(XQ));
N1=zeros(size(XQ));N2=zeros(size(XQ));
for EM=1:1000
p25=1-p15;
P1=p15*exp(-(XQ-mu15).^2/sig5^2)./(p15*exp(-(XQ-mu15).^2/sig5^2)+p25*exp(-(XQ-mu25).^2/sig5^2));
P2=p25*exp(-(XQ-mu25).^2/sig5^2)./(p15*exp(-(XQ-mu15).^2/sig5^2)+p25*exp(-(XQ-mu25).^2/sig5^2));
N1=sum(P1);N2=sum(P2);S1=sum(XQ.*P1);S2=sum(XQ.*P2);
p15=N1/n;p25=N2/n;
%mu15=S1/max(1,N1);
mu25=S2/max(1,N2);
sig5=sqrt((1/n)*(M2+mu15^2*N1+mu25^2*N2-2*(mu15*S1+mu25*S2)));
%disp([mu1 mu2 mu3 p1 1-p1-p3 p2 sig])
end
clear mu12 v12 pt p12
mumix1(j)=mu15;mumix2(j)=mu25;pmix1(j)=p15;pmix2(j)=p25;sigmix(j)=sig5;
[mu12 v12]=sort([mumix1(j) mumix2(j)]);pt=[p15 p25];p12=pt(v12);
valleymix(j)=(mumix1(j)+mumix2(j))/2-(sigmix(j)^2/(mu12(2)-mu12(1)))*log(pt(2)/pt(1));
zdistmix(j)=(mu12(2)-mu12(1))/sigmix(j);
z1=(valleymix(j)-mu12(1))/sigmix(j);z2=(-valleymix(j)+mu12(2))/sigmix(j);
z1=sign(z1)*min(abs(z1),4);z2=sign(z2)*min(abs(z2),4);
fz1=normcdf(z1)-z1*normpdf(z1);fz2=normcdf(z2)-z2*normpdf(z2);
stdestmixBOOT(j)=stdestBOOT(j,2);other(j)=0;
if zdistmix(j)>1.5 & YY(j)> valleymix(j)
stdestmixBOOT(j)=std(XQ(XQ>valleymix(j)))/fz2;
%stdestmixBOOT(j)=std(XQ(XQ>valleymix(j)));
other(j)=mean(XQ(XQ<=valleymix(j)));
end
if zdistmix(j)>1.5 & YY(j)< valleymix(j)
stdestmixBOOT(j)=std(XQ(XQ<valleymix(j)))/fz1;
%stdestmixBOOT(j)=std(XQ(XQ<valleymix(j)));
other(j)=mean(XQ(XQ>=valleymix(j)));
end
end


%RT=SHest1(1,:);for i=2:asd(2)-1,RT=[RT SHest1(i,:)];end
%mu10=mean(RT)-std(RT);mu20=mean(RT);mu30=mean(RT)+std(RT);
%p10=1/3;p30=1/3;
%sig0=std(RT)/2;
%mixture3
%clear u v pp,pp=[p1 p2 p3];[umu vmu]=sort([mu1 mu2 mu3]);upp=pp(vmu);
%valley=(mu1+mu2)/2-(sig^2/(umu(2)-umu(1)))*log(upp(2)/upp(1));

%for j=1:asd(2)-1
%    if YY(j)>valley
%        clear humpdat,humpdat=sort(SHest1(j,SHest1(j,:)>=valley));
%    end
%    if YY(j)<valley
%        clear humpdat,humpdat=sort(SHest1(j,SHest1(j,:)<=valley));
%    end
%    stdestnew(j)=std(max(humpdat(3),min(humpdat(end-2),humpdat)));
%    end

for jh=1:asd(2)-1
    ZscoreestmixBOOT(jh)=(SHAPLEYfullestBOOT(UU,jh)-mm)/stdestmixBOOT(jh);
    pvalestmixBOOT(jh)=2*(1-tcdf(abs(ZscoreestmixBOOT(jh)),asd(1)-1));
end
clear u v uu uuu
[u v]=sort(pvalestmixBOOT);
uu=u*(asd(2)-1)./(1:asd(2)-1);
uuu(1)=uu(1);
for i=2:asd(2)-1,uuu(i)=max(uuu(i-1),uu(i));end
for i=1:asd(2)-1,pvalestFDRmixBOOT(v(i))=uuu(i);end

pvalmixBOOT=[pvalestmixBOOT' pvalestFDRmixBOOT'];
SHAPandSTDmixBOOT=[SHAPLEYfullestBOOT(UU,:)' stdestmixBOOT(:)]; 
%CInew=[SHAPLEYfullest(UU,:)'-tinv(1-alpha/2,asd(1)-1)*stdestnew(:) SHAPLEYfullest(UU,:)' SHAPLEYfullest(UU,:)'+tinv(1-alpha/2,asd(1)-1)*stdestnew(:) ones(asd(2)-1,1)*mm];
CImixBOOT=[SHAPLEYfullestBOOT(UU,:)'-tinv(1-alpha/2,asd(1)-1)*stdestmixBOOT(:) SHAPLEYfullestBOOT(UU,:)' SHAPLEYfullestBOOT(UU,:)'+tinv(1-alpha/2,asd(1)-1)*stdestmixBOOT(:) ones(asd(2)-1,1)*mm];

YY=SHAPLEYfullestBOOT(UU,:);gausskervalley
CIcalibmixSHAPLBOOT=(100/(mm*(asd(2)-1)))*[calibYY'-tinv(1-alpha/2,asd(1)-1)*stdestmixBOOT(:)*factor1 calibYY' calibYY'+tinv(1-alpha/2,asd(1)-1)*stdestmixBOOT(:)*factor1 ones(asd(2)-1,1)*mm];
clear rt2 XQ2 
rt2=SHest1BOOT(1,:);for i=2:asd(2)-1,rt2=[rt2 SHest1BOOT(i,:)];end
XQ2=rt2(10:end-9);
mu10=mean(XQ2)-std(XQ2);mu20=mean(XQ2);mu30=mean(XQ2)+std(XQ2);
p10=1/3;p30=1/3;sig0=std(XQ2)/2;
mixture3
mumin=min([mu1 mu2 mu3]);factorBOOT=1/(1-mumin/aver);
calibYYBOOT=aver+factorBOOT*(YY-aver);
%disp([mu1 mu2 mu3 p1 1-p1-p3 p3 sig])
CIcalibmixBOOT=(100/(mm*(asd(2)-1)))*[calibYYBOOT'-tinv(1-alpha/2,asd(1)-1)*stdestmixBOOT(:)*factorBOOT calibYYBOOT' calibYYBOOT'+tinv(1-alpha/2,asd(1)-1)*stdestmixBOOT(:)*factorBOOT ones(asd(2)-1,1)*mm];

