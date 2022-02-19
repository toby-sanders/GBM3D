function I = LTsharpen(I,alpha,sigma)

% BM3D colored denoising of sharpened image

% Laplacian filter (high pass) is added to the original image with constant
% alpha, which enhances the edges/feature, but also brings out any noise.
% BM3D colored denoising is then used to remove then enhanced noise

% Of note: a two-step variation of the sharpening was tried, but the
% results were no better. The two step version converted the sharpening to
% an equivalent Wiener deconvolution problem and used the BM3D deblurring
% algorithm. Since this is simpler, it was kept as the best solution

% get size of image and stretch sigma to be a PSD if necessary
[m,n] = size(I);
if numel(sigma)==1
    sigma = ones(m,n)*sigma^2;
end

% sharpen the image and construct the new PSD
T = [0 -1 0; -1 4 -1; 0 -1 0];
I = I + alpha*imfilter(I,T,'replicate');
Lfilt = 1 + alpha*abs(fft2(T,m,n));
PSD = (sigma).*abs(Lfilt).^2;

% denoise with BM3D
opts.profile = 'default';
I = GBM3D_distributed(I,PSD,opts);


