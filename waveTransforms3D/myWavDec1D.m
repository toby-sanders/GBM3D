function C = myWavDec1D(Psi,Phi,levels,U,dim)

% orthogonal wavelet transform
% i.e. wavelet decomposition algorithm

if nargin<5
    dim = 3;
end

[p,q,r,s,t] = size(U);
C = zeros(p,q,r,s,t);

cnt = 0;
for i = 1:levels % loop for computing detail coefficients
    % high pass filter and downsample
    FU = fft(U,size(U,dim),dim);
    d = ifft(FU.*fft(Psi,size(U,dim),dim),size(U,dim),dim);
    C(:,:,cnt+1:cnt+size(d,dim)/2,:,:) = real(d(:,:,1:2:end,:,:));
    cnt = cnt + size(d,dim)/2;
    
    % low pass filter, downsample, and move to next level
    % U = circshift(imfilter(U,Lo_d','conv','circular'),N/2);
    U = real(ifft(FU.*fft(Phi,size(U,dim),dim),size(U,dim),dim));
    U = U(:,:,1:2:end,:,:);
end
C(:,:,cnt+1:end,:,:) = U; % last step returns approximation coefficients
