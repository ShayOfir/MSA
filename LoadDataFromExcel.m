
function  [xy, subjects, regions] = LoadDataFromExcel (xlfname, sheetname, treat_zero, varargin)
% Create a matrix for MSA from Excel file
% Structre of the excel spreadsheet:
%
%       1st Row:  coloumn names (region names + performance measure)
%       2nd to last Rows:  first subect (first cell - subject name; all other cell
%       except last - numbers between 0 to 100 corresponding to percent of
%       damage; Last column - performance measure (any number).
%       
%
%       
% Use:
%
%           [xy, subjects, regions] = LoadDataFromExcel (xlfname, sheetname, treat_zero [,minPerformance, maxPerformance])
%
% input:
%         xlfname - file name of excel file:
%         sheetname - sheet name, put 1 for the first.
%         treat_zero - What to do with regions with zero damage accross all
%         patients. 0=Don't exclude, 1=exclude
%         minPerformance - optional. if stated, it is the theoretical
%         performance measure of a completely damaged patient.
%         maxPerformance - optional. if srared it is the theoretical performance measure of the subject if neurologically-intact  
%         
%
% Output:
%          xy - matrix of m rows x n columns
%               each row correspons to a subject
%               each column from 1:end-1 correspons to a damage extent in
%               the a brain-region
%               the last column correspons to performance measure
%        
%         subjects - cell-vector of names of subjects (row names)
%         regions - cell-vector of names of regions (coloumn names)

minP = nan;
maxP = nan;
tab = readtable(xlfname, 'Sheet', sheetname,'ReadRowNames',true);
subjects = tab.Properties.RowNames';
regions = tab.Properties.VariableNames(1:end-1);
xy = table2array(tab);
x = xy(:,1:end-1);
y = xy(:,end);
if treat_zero==1  
    nonzero = x(:,sum(x)>0);
    xy = [nonzero,y];
    regions = regions(sum(x)>0);
end


if length(varargin) > 0
    minP = varargin{1};
    if length(varargin) > 1
        maxP = varargin{2};
    end
    y = xy(:,end);
    Y = y;
    if ~isnan(minP)    
        Y(y<minP) = minP;   
    end
    if ~isnan(maxP)
        Y(y>maxP) = maxP;
    end
    xy(:,end) = Y;
end




end