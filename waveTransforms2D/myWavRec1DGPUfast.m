function U = myWavRec1DGPUfast(Psi,Phi,levels,C)

% Written by Toby Sanders @Lickenbrock Tech
% 10-19-2020


% 2-dimensional orthogonal conjugate wavelet transform
% i.e. wavelet reconstruction algorithm

Cup = cell(levels+1,1);
% HPF the upsampled detail coefficients
for i = 1:levels
    Cup{i} = gpuArray(zeros(2*size(C{i,1},1),1,'single'));
    tmp = Cup{i};   
    tmp(1:2:end) = C{i};
    Cup{i} = Cup{i} + real(ifft(fft(tmp).*conj(fft(Psi,size(tmp,1)))));

end
Cup{levels+1} = C{levels+1};

% iteratively LPF and upsample
for i = levels+1:-1:2 % start at lowest levels and work up
    % upsample at current level
    tmp = gpuArray(zeros(numel(Cup{i-1}),1,'single'));
    tmp(1:2:end) = Cup{i};
    
    % filter and accumulate upsampled into next level up
    FPsi = conj(fft(Phi,numel(tmp)));  
    Cup{i-1} = Cup{i-1} + real(ifft(fftn(tmp).*FPsi)); 
end
% at end all filtered coefficients have been collected into highest level
U = Cup{1};