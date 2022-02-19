function [U,W] = denoiseAndAggregateWie(V1,V2,S,sigma,opts)

% Written by Toby Sanders @Lickenbrock Tech
% 3-2-2021


% empirical Wiener collaborative filtering of matched blocks
% inputs:
% V1 - 5D denoised volumes of matched blocks
% V2 - 5D noisy volumes of matched blocks
% S - indices of matched blocks
% sigma - noise variance
% opts - structure with options

% output: U - denoised image
useGPU = 0;

% build the volume to denoise based on matched blocks
[M,N] = size(S); % number of reference blocks
d = opts.blockSizeWie;
d2 = opts.matchSize; % window size of local search in block matching
m = M*d; n = N*d; % image size
levels = 3; % levels in wavelet transform
if size(V1,3)<16, levels = 2; end

% Haar wavelet filters
Psi = single(reshape([-sqrt(2)/2,sqrt(2)/2],1,1,2));
Phi = abs(Psi);
if gpuDeviceCount>useGPU 
    Psi = gpuArray(Psi);Phi = gpuArray(Phi);
end

% 2D DCT transform + 1D Haar transform
V1 = dct(V1,[],1,'Type',2);
V1 = dct(V1,[],2,'Type',2);
if gpuDeviceCount>useGPU , V1 = gpuArray(V1); end
V1 = myWavDec1D(Psi,Phi,levels,V1,3);

% empirical Wiener filter estimate
filt = abs(V1).^2./(abs(V1).^2 + sigma^2);
filt(1,1,15:end,:,:) = 1;

% transform noisy volume and apply Wiener filter
V2 = dct(V2,[],1,'Type',2);
V2 = dct(V2,[],2,'Type',2);
if gpuDeviceCount>useGPU, V2 = gpuArray(V2); end
V2 = myWavDec1D(Psi,Phi,levels,V2,3);
V2 = V2.*filt;

% inverse transform to obtain denoise blocks
V2 = myWavRec1D(Psi,Phi,levels,V2,3);
if gpuDeviceCount>useGPU , V2 = gather(V2); end
V2 = idct(V2,[],2,'Type',2);
V2 = idct(V2,[],1,'Type',2);
V2 = min(max(V2,0),1);

% aggregate results, looping over each block
U = zeros(m+d2,n+d2); % pad for block matches that wrap around
W = zeros(m+d2,n+d2);
% out.Z = zeros(m,n);

w = linspace(.2,.8,d);
w = sin(pi*w);
% w = kaiser(d,2);
w = w'*w;
% w(:) = 1;

for y = 1:M
    for x = 1:N
        blocks = V2(:,:,:,y,x);    
        weight = w/norm(col(filt(:,:,:,y,x)))^2;
        for i = 1:numel(S{y,x})
            % get indices and shift according to reference location
            [a,b] = ind2sub([d2,d2],S{y,x}(i));
            a = a + (y-1)*d; 
            b = b + (x-1)*d;
            
            % aggregate into denoised image and add in weights
            U(a:a+d-1,b:b+d-1) = U(a:a+d-1,b:b+d-1) + weight.*blocks(:,:,i);
            W(a:a+d-1,b:b+d-1) = W(a:a+d-1,b:b+d-1) + weight;
        end
        % out.Z((y-1)*d+1:y*d,(x-1)*d+1:x*d) = numel(S{y,x});
    end
end


% reweight image by number of times each pixel is visited
% U = U./W;
U = U(d2/2+1:m+d2/2,d2/2+1:n+d2/2); % delete the padding
W = W(d2/2+1:m+d2/2,d2/2+1:n+d2/2);