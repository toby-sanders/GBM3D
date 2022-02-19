function U = myWavRec1D(LoR,HiR,levels,C)

% orthogonal conjugate wavelet transform
% i.e. wavelet reconstruction algorithm
p = 2*size(C{1},1);

N = numel(LoR);

U = zeros(p,1); % store restored signal into U
LoR = LoR';
HiR = HiR';


Cup = cell(levels+1,1);
% HPF the upsampled detail coefficients
for i = 1:levels
    Cup{i} = zeros(2*numel(C{i}),1);
    Cup{i}(1:2:end) = C{i};
    Cup{i} = circshift(imfilter(Cup{i},HiR,'corr','circular'),-N/2+1);
end
Cup{levels+1} = C{levels+1};

% iteratively LPF and upsample
for i = 1:levels
    U = U + Cup{i}(1:p); % next level coefficients added
   % upsample coefficients and LPF
   for j = i+1:levels+1
       tmp = zeros(numel(Cup{j})*2,1);
       tmp(1:2:end) = Cup{j};
       Cup{j} =  circshift(imfilter(tmp,LoR,'corr','circular'),-N/2+1);
   end
end
   U = U + Cup{levels+1}(1:p); % next level coefficients added