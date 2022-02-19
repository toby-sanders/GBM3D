clear;
d = 128; % image dim.
sZ = 32; % size of blocks
SNR = 6; % SNR of noisy image
opts.blockSize = sZ;
opts.numMax = 16; % max number of blocks to match for each ref. block
opts.numMin = 16; % min number of blocks to match for each ref. block
opts.wname = 'bior';
opts.wnamez = 'db';
opts.order = 13;
opts.orderz = 1;
opts.levels = 3;
opts.cycleSpin = 2;
opts.matchSpins =  0;% [0 4 sZ/2-1];
opts.wiener = false;
opts.tauMode = 1;
rng(2021);

% get image and add noise
% I = imresize(im2double(rgb2gray(imread('peppers.png'))),[d,d]);
% I = im2double(imread('cameraman.tif'));
% I = im2double(imread('images/barbara.png'));
% I = im2double(rgb2gray(imread('images/IMG_0086.JPG')));
% I = I(1201:1201+d-1,1601:1601+d-1);
% I = phantom(d);
I = im2double(rgb2gray(imread('images/lena.png')));
% I = im2double(rgb2gray(imread('images/fabio.jpeg')));
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

IO = I;
[I,sigma] = add_Wnoise(I,SNR);

% our BM3D code
gpuID = gpuDevice(2);
reset(gpuID);


tic;
[U] = LTBM3D(I,sigma,opts);
toc;

tic;
[U2] = LTBM3D2(I,sigma,opts);
toc;


myrel(U,U2)