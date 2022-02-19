function [Z,M] = GBM3D_decompose(I,sigma0,levels)


[m,n,C] = size(I);
Z = zeros(m,n,C,levels+1);
Z(:,:,:,1) = I;
BMopts.profile = 'default';
if numel(sigma0)==1
    sigmas = linspace(sigma0,2*sigma0,levels);
else
    sigmas = linspace(sigma0(1),sigma0(2),levels);
end

for i = 1:levels
    for j = 1:C
        Z(:,:,j,i+1) = GBM3D_distributed(Z(:,:,j,i),sigmas(i),BMopts);
    end
end

% get the masking effect using block matching metrics
M = getMaskImage2(I(:,:,1));
M = max(M/max(M(:)),0);

