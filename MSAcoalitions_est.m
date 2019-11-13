clear XYZ XYZsize T1 T2 T3 T4 T5 T6 T7 Xdat Vdat V XX W dista1
clear VVest preSHest pre1SHest SHest Vest V1est MMpr MMpr1
b=15; BB=5;CC=0;JJJ=0;

www=size(datum); nn=www(1);
XYZ=[(datum(1:nn,1:end-1)./(ones(nn,1)*max(datum(1:nn,1:end-1)))) (datum(1:nn,end))];

% Input data XYZ of fraction damaged. One more column for performance.
% Zero-performance added for up to 3 working regions. Full-performance
% added for no damage.

XYZsize=size(XYZ);KK=XYZsize(1);m=XYZsize(2)-1;n=KK+2;
Xdat=zeros(n,m);Vdat=zeros(n,1);
V=zeros(n,1);Vest=zeros(n,1);V1est=V;
Xdat(1:KK,:)=1-XYZ(:,1:end-1);

Vdat=XYZ(:,end);
TOP=max(Vdat);
%TOP=4.5;
if III==1, outdatum=[1-XYZ(:,1:end-1) XYZ(:,end)];end
Xdat(KK+1,:)=zeros(1,m);Vdat(KK+1)=0;
Xdat(KK+2,:)=ones(1,m);Vdat(KK+2)=TOP;
% Xdat is region condition data, Vdat is subject performance data.
MMpr=fitrensemble(Xdat,Vdat);

% outdatum analysis
if III==0
XX=outdatum(qqq,1:end-1);
for j=1:KK+2
    dista1(j)=sum((Xdat(j,:)-XX).^2);dista1(j)=min(BB,dista1(j)-CC);
    W(j)=exp(-b*dista1(j));
end
    gam1=((W*Xdat)*XX')/((W*Xdat)*(W*Xdat)');
    gam2=sum(XX)/sum(W*Xdat);
    gam3=min(1,1/max(gam2*W*Xdat));
    bet=1/sum(W);
    V1est=gam3*gam2*sum(W.*Vdat');
    MMpr1=fitrensemble(Xdat,Vdat);
    
    predout(qqq)=min(TOP,gam2*sum(W.*Vdat'));
if JJJ==1,predout(qqq)=predict(MMpr1,XX);end
end


%coalitions of size m-1
IJ=1;
for i1=1:m
XX=ones(1,m);XX(i1)=0;    
for j=1:KK+2
    dista1(j)=sum((Xdat(j,:)-XX).^2);dista1(j)=min(BB,dista1(j)-CC);
    W(j)=exp(-b*dista1(j));
end
    gam1=((W*Xdat)*XX')/((W*Xdat)*(W*Xdat)');
    gam2=sum(XX)/sum(W*Xdat);
    gam3=min(1,1/max(gam2*W*Xdat));
    bet=1/sum(W);
    V1est=gam3*gam2*sum(W.*Vdat');
%V1est=min(TOP,gam2*sum(W.*Vdat'));
        
    if JJJ==1,V1est=predict(MMpr,XX);end
    T1(IJ,:)=[i1 V1est];
    IJ=IJ+1;
end

%coalitions of size m-2
IJ=1;
for i1=1:m-1,for i2=i1+1:m
XX=ones(1,m);XX(i1)=0;XX(i2)=0;    
for j=1:KK+2
    dista1(j)=sum((Xdat(j,:)-XX).^2);dista1(j)=min(BB,dista1(j)-CC);
    W(j)=exp(-b*dista1(j));
end
    gam1=((W*Xdat)*XX')/((W*Xdat)*(W*Xdat)');
    gam2=sum(XX)/sum(W*Xdat);
    gam3=min(1,1/max(gam2*W*Xdat));
    bet=1/sum(W);
    V1est=gam3*gam2*sum(W.*Vdat');
    if JJJ==1,V1est=predict(MMpr,XX);end
    T2(IJ,:)=[i1 i2 V1est];
    IJ=IJ+1;
end,end

%coalitions of size m-3
IJ=1;
for i1=1:m-2,for i2=i1+1:m-1;for i3=i2+1:m
XX=ones(1,m);XX(i1)=0;XX(i2)=0;XX(i3)=0;    
for j=1:KK+2
    dista1(j)=sum((Xdat(j,:)-XX).^2);dista1(j)=min(BB,dista1(j)-CC);
    W(j)=exp(-b*dista1(j));
end
gam1=((W*Xdat)*XX')/((W*Xdat)*(W*Xdat)');
    gam1=((W*Xdat)*XX')/((W*Xdat)*(W*Xdat)');
    gam2=sum(XX)/sum(W*Xdat);
    gam3=min(1,1/max(gam2*W*Xdat));
    bet=1/sum(W);
    V1est=gam3*gam2*sum(W.*Vdat');
    if JJJ==1,V1est=predict(MMpr,XX);end    
    T3(IJ,:)=[i1 i2 i3 V1est];
    IJ=IJ+1;
end,end,end

%coalitions of size m-4
IJ=1;
for i1=1:m-3,for i2=i1+1:m-2;for i3=i2+1:m-1;for i4=i3+1:m
XX=ones(1,m);XX(i1)=0;XX(i2)=0;XX(i3)=0;XX(i4)=0;    
for j=1:KK+2
    dista1(j)=sum((Xdat(j,:)-XX).^2);dista1(j)=min(BB,dista1(j)-CC);
    W(j)=exp(-b*dista1(j));
end
    gam1=((W*Xdat)*XX')/((W*Xdat)*(W*Xdat)');
    gam2=sum(XX)/sum(W*Xdat);
    gam3=min(1,1/max(gam2*W*Xdat));
    bet=1/sum(W);
    V1est=gam3*gam2*sum(W.*Vdat');
    if JJJ==1,V1est=predict(MMpr,XX);end
    T4(IJ,:)=[i1 i2 i3 i4 V1est];
    IJ=IJ+1;
end,end,end,end
if UU>4
%coalitions of size m-5
IJ=1;
for i1=1:m-4,for i2=i1+1:m-3,for i3=i2+1:m-2,for i4=i3+1:m-1,for i5=i4+1:m
XX=ones(1,m);XX(i1)=0;XX(i2)=0;XX(i3)=0;XX(i4)=0;XX(i5)=0;    
for j=1:KK+2
    dista1(j)=sum((Xdat(j,:)-XX).^2);dista1(j)=min(BB,dista1(j)-CC);
    W(j)=exp(-b*dista1(j));
end
    gam1=((W*Xdat)*XX')/((W*Xdat)*(W*Xdat)');
    gam2=sum(XX)/sum(W*Xdat);
    gam3=min(1,1/max(gam2*W*Xdat));
    bet=1/sum(W);
    V1est=gam3*gam2*sum(W.*Vdat');
    if JJJ==1,V1est=predict(MMpr,XX);end
    T5(IJ,:)=[i1 i2 i3 i4 i5 V1est];
    IJ=IJ+1;
end,end,end,end,end
end
if UU>5
%coalitions of size m-6
IJ=1;
for i1=1:m-5,for i2=i1+1:m-4,for i3=i2+1:m-3,for i4=i3+1:m-2,for i5=i4+1:m-1;for i6=i5+1:m
XX=ones(1,m);XX(i1)=0;XX(i2)=0;XX(i3)=0;XX(i4)=0;XX(i5)=0;XX(i6)=0;    
for j=1:KK+2
    dista1(j)=sum((Xdat(j,:)-XX).^2);dista1(j)=min(BB,dista1(j)-CC);
    W(j)=exp(-b*dista1(j));
end
    gam1=((W*Xdat)*XX')/((W*Xdat)*(W*Xdat)');
    gam2=sum(XX)/sum(W*Xdat);
    gam3=min(1,1/max(gam2*W*Xdat));
    bet=1/sum(W);
    V1est=gam3*gam2*sum(W.*Vdat');
    if JJJ==1,V1est=predict(MMpr,XX);end
    T6(IJ,:)=[i1 i2 i3 i4 i5 i6 V1est];
    IJ=IJ+1;
end,end,end,end,end,end
end
if UU>6
%coalitions of size m-7
IJ=1;
for i1=1:m-6,for i2=i1+1:m-5,for i3=i2+1:m-4,for i4=i3+1:m-3,for i5=i4+1:m-2;for i6=i5+1:m-1,for i7=i6+1:m
XX=ones(1,m);XX(i1)=0;XX(i2)=0;XX(i3)=0;XX(i4)=0;XX(i5)=0;XX(i6)=0;XX(i7)=0;    
for j=1:KK+2
    dista1(j)=sum((Xdat(j,:)-XX).^2);dista1(j)=min(BB,dista1(j)-CC);
    W(j)=exp(-b*dista1(j));
end
    gam1=((W*Xdat)*XX')/((W*Xdat)*(W*Xdat)');
    gam2=sum(XX)/sum(W*Xdat);
    gam3=min(1,1/max(gam2*W*Xdat));
    bet=1/sum(W);
    V1est=gam3*gam2*sum(W.*Vdat');
    if JJJ==1,V1est=predict(MMpr,XX);end
    T7(IJ,:)=[i1 i2 i3 i4 i5 i6 i7 V1est];
    IJ=IJ+1;
end,end,end,end,end,end,end
end
for i=1:m
clear QQQ
QQQ=T1(T1(:,1)~=i,:);
    VVest(1,i)=mean(QQQ(:,end));
    pre1SHest(1,i)=VVest(1,i);
    preSHest(1,i)=pre1SHest(1,i)-VVest(1,i)/m;
    
clear QQQ
QQQ=T2(T2(:,1)~=i & T2(:,2)~=i,:);
    VVest(2,i)=mean(QQQ(:,end));
    pre1SHest(2,i)=pre1SHest(1,i)+VVest(2,i)/2;
    preSHest(2,i)=pre1SHest(2,i)-VVest(2,i)/m;
clear QQQ
QQQ=T3(T3(:,1)~=i & T3(:,2)~=i & T3(:,3)~=i,:);
    VVest(3,i)=mean(QQQ(:,end));
    pre1SHest(3,i)=pre1SHest(2,i)+VVest(3,i)/3;
    preSHest(3,i)=pre1SHest(3,i)-VVest(3,i)/m;

clear QQQ
QQQ=T4(T4(:,1)~=i & T4(:,2)~=i & T4(:,3)~=i & T4(:,4)~=i,:);
    VVest(4,i)=mean(QQQ(:,end));
    pre1SHest(4,i)=pre1SHest(3,i)+VVest(4,i)/4;
    preSHest(4,i)=pre1SHest(4,i)-VVest(4,i)/m;
if UU>4
clear QQQ
QQQ=T5(T5(:,1)~=i & T5(:,2)~=i & T5(:,3)~=i & T5(:,4)~=i & T5(:,5)~=i,:);
    VVest(5,i)=mean(QQQ(:,end));
    pre1SHest(5,i)=pre1SHest(4,i)+VVest(5,i)/5;
    preSHest(5,i)=pre1SHest(5,i)-VVest(5,i)/m;
end
if UU>5
clear QQQ
QQQ=T6(T6(:,1)~=i & T6(:,2)~=i & T6(:,3)~=i & T6(:,4)~=i & T6(:,5)~=i  & T6(:,6)~=i,:);
    VVest(6,i)=mean(QQQ(:,end));
    pre1SHest(6,i)=pre1SHest(5,i)+VVest(6,i)/6;
    preSHest(6,i)=pre1SHest(6,i)-VVest(6,i)/m;
end
if UU>6
clear QQQ
QQQ=T7(T7(:,1)~=i & T7(:,2)~=i & T7(:,3)~=i & T7(:,4)~=i & T7(:,5)~=i  & T7(:,6)~=i & T7(:,7)~=i,:);
    VVest(7,i)=mean(QQQ(:,end));
    pre1SHest(7,i)=pre1SHest(6,i)+VVest(7,i)/7;
    preSHest(7,i)=pre1SHest(7,i)-VVest(7,i)/m;
end
end
for j=1:UU
    meanest(j)=mean(preSHest(j,:));
    SHest(j,:)=TOP/m+preSHest(j,:)-meanest(j);
end

