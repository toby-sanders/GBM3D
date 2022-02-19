clear;
d = 128; % image dim.
sZ = 32; % size of blocks
SNR = 10; % SNR of noisy image
opts.blockSize = sZ;
opts.numMax = 16; % max number of blocks to match for each ref. block
opts.numMin = 16; % min number of blocks to match for each ref. block
opts.wname = 'bior';
opts.wnamez = 'db';
opts.order = 13;
opts.orderz = 1;
opts.levels = 3;
opts.cycleSpin = 2;
opts.matchSpins =  [0 4 sZ/2-1];
opts.wiener = false;
opts.tauMode = 2;

pp = 1/(SNR+1);
rng(2021);

% get image and add noise
% I = imresize(im2double(rgb2gray(imread('peppers.png'))),[d,d]);
% I = im2double(imread('cameraman.tif'));
% I = im2double(imread('images/barbara.png'));
% I = im2double(rgb2gray(imread('images/IMG_0086.JPG')));
% I = I(1201:1201+d-1,1601:1601+d-1);
% I = phantom(d);
% I = im2double(rgb2gray(imread('images/lena.png')));
I = im2double(rgb2gray(imread('images/fabio.jpeg')));
% I = im2double(rgb2gray(imread('images/peppers2.png')));
% I = im2double(rgb2gray(imread('images/monarch.png')));
% I = im2double(rgb2gray(imread('images/tulips.png')));
% I = im2double(rgb2gray(imread('images/cat.png')));
% I = im2double(rgb2gray(imread('images/baboon.png')));
% I = ones(d)*1/2;
% I(3*d/8:5*d/8,3*d/8:5*d/8) = 1;
% I = imread('images/house.tif');
% I = im2double(I(:,:,1));
% I = im2double((imread('images/surfer.png')));

I0 = I;
% [I,sigma] = add_Wnoise(I,SNR);
[d1,d2] = size(I);
sigma = mean(abs(I(:)))/SNR;
S = randn(d1,d2);
[~,S] = sort(S(:));
S1 = S(1:round(pp*d1*d2/2));
S2 = S(round(pp*d1*d2/2)+1:round(pp*d1*d2));
I(S1) = 1;
I(S2) = 0;
% sigma = .02;

% our BM3D code
% gpuID = gpuDevice(2);
% reset(gpuID);
opts.filtType = 'median';
tic;
[U] = LTBM3D(I,sigma,opts);
toc;

opts.filtType = 'ht';
tic;
U2 = LTBM3D(I,sigma,opts);
toc;

myrel(U,U2)

figure(133);colormap(gray);
tiledlayout(2,3,'tilespacing','none');
t1 = nexttile;imagesc(I0,[0 1]);
t2 = nexttile;imagesc(I,[0 1]);
t3 = nexttile;imagesc(U,[0 1]);
t4 = nexttile;imagesc(U2,[0 1]);