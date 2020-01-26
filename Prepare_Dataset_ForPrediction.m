function xp = Prepare_Dataset_ForPrediction (x, varargin)
% Prepares raw data for prediction: normalizes damage, transform it to extent of activity and adds intact and zero patients
% 
%      xp=PrepareDatasetForPrediction(x, normalize)
%
% x - matrix of n patients X (m regions + 1 behavioral score)
% normalize - 0= Do no normalize, 1 = normalize (default)
% xp - matrix of n+2 pateints X (m regions + 1 behavioral score) 
[n, m] = size(x);
m = m - 1; 
if isempty(varargin) 
    normalize = 1;
else
    normalize = varargin{1};
end
if normalize 
    x=[x(:,1:end-1)./(ones(n,1)*max(x(:,1:end-1))), x(1:n,end)];%Normalizes damage data
end

x(isnan(x)) = 0;
x(:,1:end-1) = 1 - x(:,1:end-1);
TOP = max(x(:,end));
xp = x;
xp(n+1,:) = [zeros(1,m), 0];
xp(n+2,:) = [ones(1,m), TOP];
end