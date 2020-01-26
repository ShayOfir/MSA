

function [V1est, dista1] = ApplyPredictor (XX, Xdat, Vdat)
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
%WXdat = W * Xdat;
%gam1 = (WXdat * XX') / (WXdat * WXdat');
gam2 = sum(XX)/sum(W*Xdat);
gam3 = min(1,1/max(gam2*W*Xdat));
sum_W = sum(W); 
bet=1/sum_W;
%sumWVdat_t = sum(W.*Vdat'); 
V1est=gam3*gam2*sum(W.*Vdat');
%Vest=min(TOP,gam2*sumWVdat_t);
end