function C = myWavDec3D(wname,levels,U)

% 3-dimensional orthogonal wavelet transform
% i.e. wavelet decomposition algorithm
if strcmp(wname,'db1'), load db1Filters;
elseif strcmp(wname,'db2'), load db2Filters;
elseif strcmp(wname,'db3'), load db3Filters;
end
N = numel(Lo_d);

% build all 8 orthogonal 3D filters from the two 1D filters
Psi = zeros(N,N,N,8);
Psi(:,:,:,1) =  Lo_d'*Lo_d.*reshape(Lo_d,1,1,N);
Psi(:,:,:,2) =  Hi_d'*Lo_d.*reshape(Lo_d,1,1,N);
Psi(:,:,:,3) =  Lo_d'*Hi_d.*reshape(Lo_d,1,1,N);
Psi(:,:,:,4) =  Lo_d'*Lo_d.*reshape(Hi_d,1,1,N);
Psi(:,:,:,5) =  Hi_d'*Hi_d.*reshape(Lo_d,1,1,N);
Psi(:,:,:,6) =  Hi_d'*Lo_d.*reshape(Hi_d,1,1,N);
Psi(:,:,:,7) =  Lo_d'*Hi_d.*reshape(Hi_d,1,1,N);
Psi(:,:,:,8) =  Hi_d'*Hi_d.*reshape(Hi_d,1,1,N);

C = cell(levels+1,7); % coef. at each level
for i = 1:levels % loop for computing detail coefficients
    % high pass filter and downsample
    for j = 1:7
        C{i,j} = circshift(imfilter(U,Psi(:,:,:,j+1),'conv','circular'),[N/2,N/2,N/2]);
        C{i,j} = C{i,j}(1:2:end,1:2:end,1:2:end);
    end
    % low pass filter, downsample, and move to next level
    U = circshift(imfilter(U,Psi(:,:,:,1),'conv','circular'),[N/2,N/2,N/2]);
    U = U(1:2:end,1:2:end,1:2:end);
end
C{levels+1,1} = U; % last step returns approximation coefficients
