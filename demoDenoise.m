% demo file for image denoising with BM3D
% written by Toby Sanders @ Lickenbrock Tech.

clear;
addpath('utilities');
addpath('source');
addpath('waveTransforms2D');
addpath('waveTransforms3D');
SNR = 5; % SNR of noisy image (will add noise)
rng(2021);

% get image and add noise
% I0 = im2double(rgb2gray(imread('fabio.jpeg')));
I0 = im2double(imread('cameraman.tif'));
sigma = mean(I0(:))/SNR;
I = I0 + randn(size(I0))*sigma;
% add_Wnoise(I0,SNR);

% Denoise with GBM3D code
tic;
% default profile options, can also set to 'fast', 'superFast', or 'accuracy'
opts.profile = 'default'; 
[U] = GBM3D(I,sigma,opts);
toc;

tic; % compare with very fast variation
opts.profile = 'superFast';
U2 = GBM3D(I,sigma,opts);
toc;


%% display results
figure(11);colormap(gray);
tiledlayout(2,2,'tilespacing','compact');
t1 = nexttile;imagesc(I0,[0 1]);
title('original image');
t2 = nexttile;imagesc(I,[0,1]);
title(sprintf('noisy image: PSNR = %g',myPSNR(I0,I,1)));
t3 = nexttile;imagesc(U,[0 1]);
title(sprintf('denoised image, default profile: PSNR = %g',myPSNR(I0,U,1)));
t4 = nexttile;imagesc(U2,[0 1]);
title(sprintf('denoised image, fast profile: PSNR = %g',myPSNR(I0,U2,1)));
linkaxes([t1 t2 t3 t4]);
