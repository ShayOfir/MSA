% This program accept one damage-performance matrix datadir matrix and
% calculates from subject LOO standard errors and p-values
alpha=.1;%UU=4;
clear datum SHAPLEYfullestLOO SHAPLEYpartestLOO SHest1LOO u v uu uuu stdestLOO pvalestLOO 
clear pvalestFDRLOO CILOO ZscoreestLOO asd pvalLOO SHAPandSTDLOO CIcalibLOO YY outdatum predout
datum=datamatrix;asd=size(datamatrix);
III=1;MSAcoalitions_est,III=0;
for j=1:UU
SHAPLEYfullestLOO(j,:)=SHest(j,:);
end
clear datum II,II=(1:asd(1))';
for qqq=1:asd(1)
    datum=datamatrix(II~=qqq,:);
if 10*floor(qqq/10)==qqq, disp(qqq),end
    MSAcoalitions_est
 for j=1:UU
    SHest1LOO(:,qqq)=SHest(j,:)';
 end
    SHAPLEYpartestLOO=mean(SHest1LOO')';
end
aa=corrcoef([outdatum(:,end) predout']);rho=aa(1,2);
for jh=1:asd(2)-1
    [u v]=sort(SHest1LOO(jh,:));
    uu=max(u(3),min(u(end-2),u));
    stdestLOO(jh,:)=[std(u) std(uu)];
end
for jh=1:asd(2)-1,
    [u v]=sort(SHest1LOO(jh,:));
    uu=max(u(3),min(u(end-2),u));
    stdestLOO(jh,:)=[std(u) std(uu)];
end
clear Z u v uu
Z=randn(asd(1),100000);
for is=1:100000
    [u v]=sort(Z(:,is));
    uu=max(u(3),min(u,u(end-2)));
    dfg(is,:)=[std(u) std(uu)];
end
var=mean(dfg.^2);
stdestLOO=stdestLOO*sqrt(asd(1)-1)/sqrt(var(2));
mm=mean(SHAPLEYfullestLOO(UU,:));
for jh=1:asd(2)-1
    ZscoreestLOO(jh)=(SHAPLEYfullestLOO(UU,jh)-mm)/stdestLOO(jh,2);
    pvalestLOO(jh)=2*(1-tcdf(abs(ZscoreestLOO(jh)),asd(1)-1));
end
clear u v uu uuu
[u v]=sort(pvalestLOO);
uu=u*(asd(2)-1)./(1:asd(2)-1);
uuu(1)=uu(1);
for i=2:asd(2)-1,uuu(i)=max(uuu(i-1),uu(i));end
for i=1:asd(2)-1,pvalestFDRLOO(v(i))=uuu(i);end
pvalLOO=[pvalestLOO' pvalestFDRLOO'];
SHAPandSTDLOO=[SHAPLEYfullestLOO(UU,:)' stdestLOO(:,2)]; 
CILOO=[SHAPLEYfullestLOO(UU,:)'-tinv(1-alpha/2,asd(1)-1)*stdestLOO(:,2) SHAPLEYfullestLOO(UU,:)' SHAPLEYfullestLOO(UU,:)'+tinv(1-alpha/2,asd(1)-1)*stdestLOO(:,2) ones(asd(2)-1,1)*mm];
YY=SHAPLEYfullestLOO(UU,:);gausskervalley
CIcalibLOO=(100/(mm*(asd(2)-1)))*[calibYY'-tinv(1-alpha/2,asd(1)-1)*stdestLOO(:,2)*factor1 calibYY' calibYY'+tinv(1-alpha/2,asd(1)-1)*stdestLOO(:,2)*factor1];
