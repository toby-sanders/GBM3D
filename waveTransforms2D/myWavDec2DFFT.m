function C = myWavDec2DFFT(Psi,levels,U)

% Written by Toby Sanders @Lickenbrock Tech
% 10-19-2020


% 2-dimensional orthogonal wavelet transform
% i.e. wavelet decomposition algorithm

C = cell(levels+1,3); % coef. at each level
for i = 1:levels % loop for computing detail coefficients
    % high pass filter and downsample
    FU = fftn(U);
    for j = 1:3
        C{i,j} = ifftn(FU.*fftn(Psi(:,:,j+1),[size(U,1),size(U,2)]));
        C{i,j} = real(C{i,j}(1:2:end,1:2:end));
    end
    % low pass filter, downsample, and move to next level
    % U = circshift(imfilter(U,Psi(:,:,:,1),'conv','circular'),[Ny/2,Nx/2,Nz/2]);
    U = real(ifftn(FU.*fftn(Psi(:,:,1),[size(U,1),size(U,2)])));
    U = U(1:2:end,1:2:end);
end
C{levels+1,1} = U; % last step returns approximation coefficients
