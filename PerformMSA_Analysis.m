% MSA FINAL CODE
function [SV, Calib, coal, d, Bset, Lset]=PerformMSA_Analysis (xy, pdepth, nBS, alpha, varargin)
%Performs MSA analysis on a dataset
%
%  PerformMSA_Analysis (prefix,pdepth,nBS,alpha,normalize)
%
% prefix - output filename with full path
% pdepth - bounded perturbation depth (1.. no. of regions); integer
% nBS -  Leave-One-Out / Boostrap:number of boostraps. 
%       > 0 = Boostrap. Than nBS is the number of bootsraps
%       = 0 = Compute SV only
%       -1 = Leave-One-Out
%
% alpha - type I error level. Default is 0.05
% 

    if isempty(varargin) 
        normalize = 1;
    else
        normalize = varargin{1} ;
    end

    [SV, coal, d, Calib] = Compute_ShapleyVector_Bound (xy, pdepth, normalize);
    Bset = cell(1,pdepth);
    Lset = cell(1,pdepth);
    if nBS > 0
        for p=pdepth
            Bset{p} = Compute_Bootstrap(xy, SV, p, nBS, alpha, normalize);
        end
    end
    if nBS == -1
        for p=pdepth           
            Lset{p} = Compute_LOO(xy, SV, p, alpha, normalize);
        end        
    end
    %save ([prefix '_' nm{ds} '.mat'],'xy','SV','Calib','coal','d','Bset','nBS','alpha','pdepth','Lset');


end


 









