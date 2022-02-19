function U = myWavRec3DGPUfast(Psi,levels,C)

% Written by Toby Sanders @Lickenbrock Tech
% 10-19-2020


% 3-dimensional orthogonal conjugate wavelet transform
% i.e. wavelet reconstruction algorithm

Cup = cell(levels+1,1);
% HPF the upsampled detail coefficients
for i = 1:levels
    Cup{i} = (zeros(2*size(C{i,1},1),2*size(C{i,1},2),2*size(C{i,1},3),'single'));
    if gpuDeviceCount>0
       Cup{i} = gpuArray(Cup{i}); 
    end
    tmp = Cup{i};
    for j = 1:7        
        tmp(1:2:end,1:2:end,1:2:end) = C{i,j};
        Cup{i} = Cup{i} + real(ifftn(fftn(tmp).*conj(fftn(Psi(:,:,:,j+1),...
            [size(tmp,1),size(tmp,2),size(tmp,3)]))));
    end
end
Cup{levels+1} = C{levels+1};

% iteratively LPF and upsample
for i = levels+1:-1:2 % start at lowest levels and work up
    % upsample at current level
    tmp = (zeros(size(Cup{i-1})));
    if gpuDeviceCount>0
        tmp = gpuArray(tmp);
    end
    tmp(1:2:end,1:2:end,1:2:end) = Cup{i};
    
    % filter and accumulate upsampled into next level up
    FPsi = conj(fftn(Psi(:,:,:,1),[size(tmp,1),size(tmp,2),size(tmp,3)]));   
    Cup{i-1} = Cup{i-1} + real(ifftn(fftn(tmp).*FPsi)); 
end
% at end all filtered coefficients have been collected into highest level
U = Cup{1};