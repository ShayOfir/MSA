
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