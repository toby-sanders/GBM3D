function U = myWavRec(wname,levels,C)

% 2-dimensional orthogonal conjugate wavelet transform
% i.e. wavelet reconstruction algorithm
[p,q] = size(C{1});
if strcmp(wname,'db1'), load db1Filters;
elseif strcmp(wname,'db2'), load db2Filters;
elseif strcmp(wname,'db3'), load db3Filters;
end
N = numel(Lo_d);
U = zeros(2*p,2*q); % store restored signal into U

% build all 4 orthogonal 2D filters from the two 1D filters
Psi = zeros(N,N,4);
Psi(:,:,1) = Lo_d'*Lo_d;
Psi(:,:,2) = Hi_d'*Hi_d;
Psi(:,:,3) = Hi_d'*Lo_d;
Psi(:,:,4) = Lo_d'*Hi_d;


Cup = cell(levels+1,3);
% HPF the upsampled detail coefficients
for i = 1:levels
    for j = 1:3
        Cup{i,j} = zeros(2*size(C{i,j},1),2*size(C{i,j},2));
        Cup{i,j}(1:2:end,1:2:end) = C{i,j};
        Cup{i,j} = circshift(imfilter(Cup{i,j},Psi(:,:,j+1),'corr','circular'),[-N/2+1,-N/2+1]);
    end
end
Cup{levels+1,1} = C{levels+1,1};

% iteratively LPF and upsample
for i = 1:levels
    U = U + Cup{i,1} + Cup{i,2} + Cup{i,3}; % next level coefficients added
    % upsample coefficients and LPF
    for k = i+1:levels+1
        for j = 1:3
            tmp = zeros(2*size(Cup{k,j},1),2*size(Cup{k,j},2));
            tmp(1:2:end,1:2:end) = Cup{k,j};
            Cup{k,j} =  circshift(imfilter(tmp,Psi(:,:,1),'corr','circular'),[-N/2+1,-N/2+1]);
        end
    end
end
U = U+Cup{levels+1,1};
