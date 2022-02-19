clear;
d = 512; % image dim.
sZ = 16; % size of blocks
SNR = 5; % SNR of noisy image

opts.numMax = 16; % max number of blocks to match for each ref. block
opts.numMin = 16; % min number of blocks to match for each ref. block
opts.wname = 'bior';
opts.wnamez = 'db';
opts.order = 15;
opts.orderz = 1;
opts.levels = 3;
opts.cycleSpin = 2;
opts.matchSpins =  [0 4 6];
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
% I = im2double(imread('images/barbara.png'));
% I = im2double(rgb2gray(imread([path,'IMG_0086.JPG'])));
% I = I(1201:1201+d-1,1601:1601+d-1);
% I = phantom(d);
% I = im2double(rgb2gray(imread([path,'lena.png'])));
% I = im2double(rgb2gray(imread([path,'fabio.jpeg'])));
% I = im2double(rgb2gray(imread([path,'peppers2.png'])));
% I = im2double(rgb2gray(imread([path,'monarch.png'])));
% I = im2double(rgb2gray(imread([path,'tulips.png'])));
I = im2double(rgb2gray(imread([path,'cat.png'])));
% I = im2double(rgb2gray(imread([path,'baboon.png'])));
% I = ones(d)*1/2;
% I(3*d/8:5*d/8,3*d/8:5*d/8) = 1;
% I = imread('images/house.tif');
% I = im2double(I(:,:,1));
% I = im2double((imread([path,'surfer.png'])));

IO = I;
[I,sigma] = add_Wnoise(I,SNR);

% our BM3D code
gpuID = gpuDevice(2);
reset(gpuID);


opts.matchSpins =  [0 4 6];
opts.blockSize = 16;
opts.blockSizeWie = 8;
opts.matchSize = 32;


tic;
[recBM1,out] = LTBM3D_distributed(I,sigma,opts);
tt1 = toc;


tic;
[recBM1,out] = LTBM3D_distributed(I,sigma,opts);
tt1 = toc;

opts.matchSpins =  [0 4];
opts.blockSize = 16;
opts.blockSizeWie = 8;
opts.matchSize = 32;

tic;
[recBM2,~] = LTBM3D_distributed(I,sigma,opts);
tt2 = toc;

opts.matchSpins =  [0 4 6];
opts.blockSize = 32;
opts.blockSizeWie = 16;
opts.matchSize = 32;

tic;
[recBM3,out] = LTBM3D_distributed(I,sigma,opts);
tt3 = toc;

opts.matchSpins =  [0 4];
opts.blockSize = 32;
opts.blockSizeWie = 8;
opts.matchSize = 32;
tic;
[recBM4,out] = LTBM3D_distributed(I,sigma,opts);
tt4 = toc;

%%
psnr1 = myPSNR(IO,recBM1,1);
psnr2 = myPSNR(IO,recBM2,1);
psnr3 = myPSNR(IO,recBM3,1);
psnr4 = myPSNR(IO,recBM4,1);



figure(1789);
subplot(1,2,1);
plot([tt1 tt2 tt3 tt4],'-o');
title('run time');xlabel('case number');
subplot(1,2,2);
plot([psnr1 psnr2 psnr3 psnr4],'-o');
xlabel('case number');
title('PSNR');









