function U = myWavRec3D_preallocated(Psi,levels,C)

% Written by Toby Sanders @Lickenbrock Tech
% 10-19-2020


% 3-dimensional orthogonal conjugate wavelet transform
% i.e. wavelet reconstruction algorithm
[p,q,r] = size(C{1});
[Ny,Nx,Nz,~] = size(Psi);
U = zeros(2*p,2*q,2*r); % store restored signal into U

Cup = cell(levels+1,7);
% HPF the upsampled detail coefficients
for i = 1:levels
    for j = 1:7
        Cup{i,j} = zeros(2*size(C{i,j},1),2*size(C{i,j},2),2*size(C{i,j},3));
        Cup{i,j}(1:2:end,1:2:end,1:2:end) = C{i,j};
        Cup{i,j} = circshift(imfilter(Cup{i,j},Psi(:,:,:,j+1),'corr','circular'),[-Ny/2+1,-Nx/2+1,-Nz/2+1]);
    end
end
Cup{levels+1,1} = C{levels+1,1};

% iteratively LPF and upsample
for i = 1:levels
    for j = 1:7
        U = U + Cup{i,j}; % next level coefficients added
    end
    % upsample coefficients and LPF
    for k = i+1:levels+1
        for j = 1:7
            tmp = zeros(2*size(Cup{k,j},1),2*size(Cup{k,j},2),2*size(Cup{k,j},3));
            tmp(1:2:end,1:2:end,1:2:end) = Cup{k,j};
            Cup{k,j} =  circshift(imfilter(tmp,Psi(:,:,:,1),'corr','circular'),[-Ny/2+1,-Nx/2+1,-Nz/2+1]);
        end
    end
end
U = U+Cup{levels+1,1};