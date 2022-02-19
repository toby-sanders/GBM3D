function U = myWavRec1D(Psi,Phi,levels,C,dim)

% orthogonal conjugate wavelet transform
% i.e. wavelet reconstruction algorithm
if nargin<5, dim = 5; end
[p,q,r,s,t] = size(C);
Cup = cell(levels+1,1);
cnt = 0;
% HPF the upsampled detail coefficients
for i = 1:levels
    Cup{i} = zeros(p,q,r/2^(i-1),s,t);
    Cup{i}(:,:,1:2:end,:,:) = C(:,:,cnt+1:cnt+r/2^i,:,:);
    cnt = cnt + r/2^i;
    FPsi = conj(fft(Psi,size(Cup{i},dim),dim));
    Cup{i} = real(ifft(fft(Cup{i},size(Cup{i},dim),dim).*FPsi,size(Cup{i},dim),dim));
end
Cup{levels+1} = C(:,:,cnt+1:end,:,:);

% iteratively LPF and upsample
for i = levels+1:-1:2
    tmp = zeros(p,q,size(Cup{i-1},dim),s,t,'single');
    tmp(:,:,1:2:end,:,:) = Cup{i};
    FPsi = conj(fft(Phi,size(tmp,dim),dim));
    Cup{i-1} = Cup{i-1} + real(ifft(fft(tmp,size(tmp,dim),dim).*FPsi,size(tmp,dim),dim));
end
U = Cup{1}; % next level coefficients added