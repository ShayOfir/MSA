function Bset = Compute_LOO (datum, SV, pdepth ,alpha, varargin)
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
UU = pdepth;
Nbig = 100000;
datamatrix = datum;
[n,m]=size(datamatrix);m=m-1;  
SHAPLEYfullest=SV; 
SHest = zeros(m,n);
Bset.stdest  = zeros(m,2);
II=(1:n)'; 
Bset.save_patients = zeros(n);
for qqq=1:n
    disp(['Computing LOO-Shapley...',int2str(qqq),'/',int2str(n)]);           
    datum=datamatrix(II~=qqq,:);
    SHest = Compute_ShapleyVector_Bound (datum, pdepth, normalize);
    SHest1LOO(:,qqq) = SHest(UU,:)';    
end
Bset.LOO = SHest1LOO;
SHAPLEYpartestLOO=mean(SHest1LOO')';
%
for jh=1:m
    [u, v]=sort(SHest1LOO(jh,:));
    uu=max(u(3),min(u(end-2),u));
    Bset.stdest(jh,:)=[std(u), std(uu)];
end
% rng(3); %<---DEBUG
Z=randn(n,Nbig);
dfg = zeros(Nbig,2);
for is=1:Nbig
    [u v]=sort(Z(:,is));
    uu=max(u(3),min(u,u(end-2)));
    dfg(is,:)=[std(u) std(uu)];
end

var=mean(dfg.^2);
Bset.stdest=Bset.stdest*sqrt(n-1)/sqrt(var(2));
mm=mean(SHAPLEYfullest(UU,:));
for jh=1:m
    Bset.Zscoreest(jh)=(SHAPLEYfullest(UU,jh)-mm)/Bset.stdest(jh,2);
    Bset.pvalest(jh)=2*(1-tcdf(abs(Bset.Zscoreest(jh)),n-1));
end
[u v]=sort(Bset.pvalest);
uu=u*m./(1:m);
uuu(1)=uu(1);
for i=2:m,uuu(i)=max(uuu(i-1),uu(i));end
for i=1:m,Bset.pvalestFDR(v(i))=uuu(i);end
Bset.pval=[Bset.pvalest' Bset.pvalestFDR'];

Bset.SHAPandSTDLOO=[SHAPLEYfullest(UU,:)' Bset.stdest(:,2)]; 
Bset.CI=[SHAPLEYfullest(UU,:)'-tinv(1-alpha/2,n-1)*Bset.stdest(:,2), SHAPLEYfullest(UU,:)', SHAPLEYfullest(UU,:)'+tinv(1-alpha/2,n-1)*Bset.stdest(:,2), ones(m,1)*mm];
YY=SHAPLEYfullest(UU,:);
[calibYY, aver, mode1, factor1] = CalibrateShapleyVector(YY);
Bset.CIcalib=(100/(mm*m))*[calibYY'-tinv(1-alpha/2,n-1)*Bset.stdest(:,2)*factor1, calibYY' calibYY'+tinv(1-alpha/2,n-1)*Bset.stdest(:,2)*factor1];


end