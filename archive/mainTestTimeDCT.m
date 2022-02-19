clear;
d = 128; % image dim.
SNR = 5; % SNR of noisy image
opts.numMax = 16; % max number of blocks to match for each ref. block
opts.numMin = 16; % min number of blocks to match for each ref. block
opts.wname = 'bior';
opts.wnamez = 'db';
opts.order = 15;
opts.orderz = 1;
opts.levels = 3;
opts.cycleSpin = 2;
opts.matchSpins =  [0 4];
opts.tauMode = 2;
opts.blockSize = 16;
opts.blockSizeWie = 8;
opts.matchSize = 32;
opts.Wiener = true;
path = 'C:\Users\toby.sanders\Dropbox\archives\data\testImages\';
rng(2021);

% get image and add noise
% I = imresize(im2double(rgb2gray(imread('peppers.png'))),[d,d]);
% I = im2double(imread('cameraman.tif'));
% I = im2double(imread([path,'barbara.png']));
% I = im2double(rgb2gray(imread('images/IMG_0086.JPG')));
% I = I(1201:1201+d-1,1601:1601+d-1);
% I = phantom(d);
I = im2double(rgb2gray(imread([path,'lena.png'])));
% I = im2double(rgb2gray(imread([path,'fabio.jpeg'])));
% I = im2double(rgb2gray(imread([path,'peppers2.png'])));
% I = im2double(rgb2gray(imread([path,'monarch.png'])));
% I = im2double(rgb2gray(imread([path,'tulips.png'])));
% I = im2double(rgb2gray(imread([path,'cat.png'])));
% I = im2double(rgb2gray(imread([path,'baboon.png'])));
% I = ones(d)*1/2;
% I(3*d/8:5*d/8,3*d/8:5*d/8) = 1;
% I = imread('images/house.tif');
% I = im2double(I(:,:,1));
% I = im2double((imread([path,'surfer.png'])));

IO = I;
[I,sigma] = add_Wnoise(I,SNR);

% our BM3D code
gpuID = gpuDevice(1);
reset(gpuID);
tic;
opts.Wiener = false;
[U] = LTBM3D(I,sigma,opts);
toc;


%% display
opts.Wiener = true;
tic;
[recBM4,out] = LTBM3D(I,sigma,opts);
toc;

%%

tic;
[recBM5,out] = LTBM3D(I,sigma,opts);
toc;

%%

myrel(recBM4,IO)
myrel(recBM5,IO)