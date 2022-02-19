% demo file for image denoising with BM3D
% written by Toby Sanders @ Lickenbrock Tech.

clear;
% addpath('utilities');
SNR = 4; % SNR of noisy image (will add noise)
% path = '/Users/tobysanders/Dropbox/archives/data/testImages/';
path = 'C:\Users\toby.sanders\Dropbox\archives\data\testImages\';
rng(2021);
% get image and add noise
% I0 = im2double(imread('cameraman.tif'));
I0 = im2double(rgb2gray(imread([path,'lena.png'])));
[I,sigma] = add_Wnoise(I0,SNR);

tic;
U = NLmeans(I,sigma);
toc;
tic;
U3 = NLmeansPatch(I,sigma);
toc;
tic;
U2 = imnlmfilt(I);
toc;
tic;
U4 = GBM3D(I,sigma);
toc;
%%
figure(122);tiledlayout(2,3,'tilespacing','none');
t1 = nexttile;imagesc(I,[0 1]);colormap(gray);
t2 = nexttile;imagesc(U,[0 1]);title(sprintf('my NL, PSNR = %g',myPSNR(I0,U,1)));
t3 = nexttile;imagesc(U3,[0 1]);title(sprintf('my NL patch, PSNR = %g',myPSNR(I0,U3,1)));
t4 = nexttile;imagesc(U2,[0 1]);title(sprintf('matlab NL, PSNR = %g',myPSNR(I0,U2,1)));
t5 = nexttile;imagesc(U4,[0 1]);title(sprintf('GBM3D, PSNR = %g',myPSNR(I0,U4,1)));
linkaxes([t1 t2 t3 t4 t5]);