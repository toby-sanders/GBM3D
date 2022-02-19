function C = myWavDec3D_preallocated(Psi,levels,U)

% Written by Toby Sanders @Lickenbrock Tech
% 10-19-2020


% 3-dimensional orthogonal wavelet transform
% i.e. wavelet decomposition algorithm
[Ny,Nx,Nz,~] = size(Psi);

C = cell(levels+1,7); % coef. at each level
for i = 1:levels % loop for computing detail coefficients
    % high pass filter and downsample
    for j = 1:7
        C{i,j} = circshift(imfilter(U,Psi(:,:,:,j+1),'conv','circular'),[Ny/2,Nx/2,Nz/2]);
        C{i,j} = C{i,j}(1:2:end,1:2:end,1:2:end);
    end
    % low pass filter, downsample, and move to next level
    U = circshift(imfilter(U,Psi(:,:,:,1),'conv','circular'),[Ny/2,Nx/2,Nz/2]);
    U = U(1:2:end,1:2:end,1:2:end);
end
C{levels+1,1} = U; % last step returns approximation coefficients
