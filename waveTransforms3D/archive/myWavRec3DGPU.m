function U = myWavRec3DGPU(Psi,levels,C)

% Written by Toby Sanders @Lickenbrock Tech
% 10-19-2020


% 3-dimensional orthogonal conjugate wavelet transform
% i.e. wavelet reconstruction algorithm
[p,q,r] = size(C{1});
U = gpuArray(zeros(2*p,2*q,2*r,'single')); % store restored signal into U

Cup = cell(levels+1,7);
% HPF the upsampled detail coefficients
for i = 1:levels
    for j = 1:7
        Cup{i,j} = gpuArray(zeros(2*size(C{i,j},1),2*size(C{i,j},2),2*size(C{i,j},3),'single'));
        Cup{i,j}(1:2:end,1:2:end,1:2:end) = C{i,j};
        Cup{i,j} = real(ifftn(fftn(Cup{i,j}).*conj(fftn(Psi(:,:,:,j+1),...
            [size(Cup{i,j},1),size(Cup{i,j},2),size(Cup{i,j},3)]))));
    end
end
Cup{levels+1,1} = C{levels+1,1};

% iteratively LPF and upsample
for i = 1:levels
    for j = 1:7
        U = U + Cup{i,j}; % next level coefficients added
        Cup{i,j} = 0;
    end
    % upsample coefficients and LPF
    for k = i+1:levels+1
        tmp = gpuArray(zeros(2*size(Cup{k,1},1),2*size(Cup{k,1},2),2*size(Cup{k,1},3),'single'));
        FPsi = conj(fftn(Psi(:,:,:,1),[size(tmp,1),size(tmp,2),size(tmp,3)]));
        for j = 1:7           
            tmp(1:2:end,1:2:end,1:2:end) = Cup{k,j};
            Cup{k,j} = real(ifftn(fftn(tmp).*FPsi));
            if k==levels+1, break; end
        end
    end
end
U = U+Cup{levels+1,1};