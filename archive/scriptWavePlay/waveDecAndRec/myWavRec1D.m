function U = myWavRec1D(wname,levels,C)

% orthogonal conjugate wavelet transform
% i.e. wavelet reconstruction algorithm
p = 2*size(C{1},1);
if strcmp(wname,'db1'), load db1Filters;
elseif strcmp(wname,'db2'), load db2Filters;
elseif strcmp(wname,'db3'), load db3Filters;
end
N = numel(Lo_d);

U = zeros(p,1); % store restored signal into U
Lo_d = Lo_d';
Hi_d = Hi_d';


Cup = cell(levels+1,1);
% HPF the upsampled detail coefficients
for i = 1:levels
    Cup{i} = zeros(2*numel(C{i}),1);
    Cup{i}(1:2:end) = C{i};
    Cup{i} = circshift(imfilter(Cup{i},Hi_d,'corr','circular'),-N/2+1);
end
Cup{levels+1} = C{levels+1};

% iteratively LPF and upsample
for i = 1:levels
    U = U + Cup{i}(1:p); % next level coefficients added
   % upsample coefficients and LPF
   for j = i+1:levels+1
       tmp = zeros(numel(Cup{j})*2,1);
       tmp(1:2:end) = Cup{j};
       Cup{j} =  circshift(imfilter(tmp,Lo_d,'corr','circular'),-N/2+1);
   end
end
   U = U + Cup{levels+1}(1:p); % next level coefficients added