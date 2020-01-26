
function [pxy, parcellation, RegionNames] = ParcelateRegions (xy, ParcTab, varargin)
%Parcelate brain regions
%
% Use:
%
%       [pxy, parcellation] = ParcelateRegions (xy, ParcTab, [DontParcellate, ParcellateToRoB])
%
% Input:
%          xy - matrix with numbers between 0 to 100 denoting damage
%          percent (double). Each row corresponds to a subject, each column
%          to a region, except of the last which correspons to the
%          performance measure
% 
%         ParcTab - a table (generally read from excel with the following
%         columns:
%            - RegionCode - number between 1 and maximal number of regions.
%            All regions in xy must be represented!
%            - SuperCode - number above zero corresponding to a super-roi,
%            which may consist one or more regions. 
%            
%            For example:
%            RegionCode     SuperCode
%            ----------     ---------
%                 1             1
%                 2             1
%                 3             2
%                 4             2
%                 5             3
%
%           5 regions are parcellated to 3 super-regions. First consists of
%           1 and 2, second consists of 3 and 4, third consists of 5.
%
%           SuperCode of 0 means - do not take.
%
%           DontParcellate - optional. vector of roi numbers that are
%           excluded from their super-rois. For the above exxample
%           DontParcellate = [1] will exclude region 1 from the second
%           super-region and will give it a new super-region (#4) with a
%           single region.
%         
%           ParcellateToRoB - optional. vector of roi numbers that are
%           excoluded from their super-rois and are assembeld in on
%           added "Rest of Brain" super-region. 
%

% X = xy(:,1:end-1);
% y = xy(:,end);
% take = sum(X)>0;
% x = X(:,take);
% ParcTab = ParcTab(take,:);

x = xy(:,1:end-1);
y = xy(:,end);

L = size(x,2);
% if size(ParcTab,1) ~= L
%     error('Size of parcellation must have the same number of regions ''xy''');
% end
DontCollapse = [];
CollapseToRest = [];

if length(varargin)>0
    DontCollapse=varargin{1};
    if length(varargin)>1
        CollapseToRest=varargin{2};
    end
end

nSuper = max(ParcTab.SuperCode);
collapse = zeros(L,nSuper);
regions = cell(1,nSuper);
for roi=1:L    
    supercode = ParcTab.SuperCode(roi);  
    if supercode > 0 
        if sum(collapse(:,supercode)) == 0
            regions{supercode} = ParcTab.SuperName{roi};
        end
        collapse(roi,supercode) = 1;                    
    end
end
% dbg = 0;
for k=1:length(DontCollapse)
    roi = DontCollapse(k);
    super = collapse(roi,:) ~= 0;
    if sum(collapse(:,super)) > 1 %roi is not single in super-roi, exclusion for collapse is requested
       %Create a new super-roi
       nSuper = nSuper + 1;
       newSuper = nSuper;
       %Cancel current allocation
       collapse(roi,super) = 0;
       %Add the new allocation
       collapse(roi,newSuper) = 1;       
    end
end

RoB = false;
for k=1:length(CollapseToRest)
    roi = CollapseToRest(k);
    super = collapse(roi,:) ~= 0;
   %Create a new super-roi: RoB
   if ~RoB
    nSuper = nSuper + 1;        
    newSuper = nSuper;
    RoB = true;
   else
       newSuper = nSuper;
   end
   %Cancel current allocation
   collapse(roi,super) = 0;
   %Add the new allocation
   collapse(roi,newSuper) = 1;         
end
% x = xy(:,1:end-1);
voxels = repmat(ParcTab.SizeInVoxels,1,size(x,1))';
x_voxels = x .* voxels / 100;
pxy = zeros(size(xy,1),nSuper+1);
pxy(:,end) = y;

for sup=1:nSuper
    dmg_voxels = sum(x_voxels(:,collapse(:,sup)==1),2); %no. of damaged voxels in the super-roi 'sup'
    dmg_denominator =  sum(ParcTab.SizeInVoxels(collapse(:,sup)==1)); %no. of of voxels in the super-roi 'sup
    pxy(:,sup) = 100*dmg_voxels / dmg_denominator;    
end

take = sum(pxy(:,1:end-1)>0)>0;
x = pxy(:,take);
pxy = [x, y];
nxt = 1;
for k=1:length(take)
    if take(k)==1
        RegionNames{nxt} = regions{k};
        nxt = nxt + 1;
    end
end
collapse = collapse(:,take);


parcellation = collapse;
end