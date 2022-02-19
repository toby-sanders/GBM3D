function C = myWavDec1DFFT(Psi,Phi,levels,U)

% Written by Toby Sanders @Lickenbrock Tech
% 10-19-2020


% 2-dimensional orthogonal wavelet transform
% i.e. wavelet decomposition algorithm

C = cell(levels+1,1); % coef. at each level
for i = 1:levels % loop for computing detail coefficients
    % high pass filter and downsample
    FU = fft(U);
    C{i} = ifft(FU.*fft(Psi,size(U,1)));
    C{i} = real(C{i}(1:2:end));
    U = real(ifft(FU.*fft(Phi,size(U,1))));
    U = U(1:2:end);
end
C{levels+1,1} = U; % last step returns approximation coefficients
