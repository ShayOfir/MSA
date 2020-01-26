

function [SV, SaveCoal, Dist, Calib] = Compute_ShapleyVector_Bound(datum, pdepth, varargin)
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

if isempty(varargin)
    normalize = 1;
else
    normalize = varargin{1};
end

datum=Prepare_Dataset_ForPrediction(datum, normalize);
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
% if pdepth > 7
%     MAX_DMG_REGIONS = 7;
% else
%     MAX_DMG_REGIONS = pdepth;
% end     
MAX_DMG_REGIONS = pdepth;
 Vest = cell(1,MAX_DMG_REGIONS);  
 alldist = cell(1,MAX_DMG_REGIONS);
 coalitions = cell(1,MAX_DMG_REGIONS);
 for nR = 1:MAX_DMG_REGIONS
     %disp(['Calculating Depth ' num2str(nR)]);
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
 pre1SHest = zeros(MAX_DMG_REGIONS,m);
 preSHest = zeros(MAX_DMG_REGIONS,m);
 for region=1:m
     for nR=1:MAX_DMG_REGIONS
         %find all coalitions in which region m is present
         indices = find(sum(coalitions{nR}~=region,2)==nR);
         %Compute mean of(Tij/k-1)          
         VVest(nR,region) = mean(Vest{nR}(indices));
         % Compute Sum[Tij/1 + Tij/2 + ... Tij/k]
         if nR == 1 
            pre1SHest(nR,region) = VVest(nR,region);
            preSHest(nR,region) = pre1SHest(nR,region) - VVest(nR,region)/m;
         else  
            pre1SHest(nR,region) = pre1SHest(nR-1,region) + VVest(nR,region)/nR;
            preSHest(nR,region) = pre1SHest(nR,region)- VVest(nR,region)/m;
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
 for j=1:size(SV,1)
    [Calib.SV(j,:), Calib.aver(j), Calib.mode1(j), Calib.factor1(j)] = CalibrateShapleyVector(SV(j,:));
 end
 
 SaveCoal = coalitions;
 Dist = alldist;   
end %elementary_MSA

