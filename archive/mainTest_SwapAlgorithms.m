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
% I = im2double(imread([path,'barbara.png']));
% I = im2double(rgb2gray(imread('images/IMG_0086.JPG')));
% I = I(1201:1201+d-1,1601:1601+d-1);
% I = phantom(d);
% I = im2double(rgb2gray(imread([path,'lena.png'])));
% I = im2double(rgb2gray(imread([path,'fabio.jpeg'])));
I = im2double(rgb2gray(imread([path,'peppers2.png'])));
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

opts.Wiener = false;
recG1 = LTBM3D(I,sigma,opts);

% denoise for baseline comparisons
profile = 'np';
recB1 = BM3D(I,sigma,profile,BM3DProfile.HARD_THRESHOLDING);

profile = 'refilter';
recGB = BM3D(I,sigma,profile,recG1);

recBG = LTBM3D_Wie(I,recB1,sigma,opts);

opts.Wiener = true;
recG2 = LTBM3D(I,sigma,opts);

recB2 = BM3D(I,sigma);


%%
figure(334);
tiledlayout(3,3,'tilespacing','none');colormap(gray);
% t1 = nexttile;imagesc(IO,[0 1]);title('original');
t2 = nexttile;imagesc(I,[0 1]);title(sprintf('noisy: PSNR = %g',myPSNR(IO,I,1)));
t3 = nexttile;imagesc(recG1,[0 1]);title(sprintf('GBM3D1: PSNR = %g',myPSNR(IO,recG1,1)));
t4 = nexttile;imagesc(recB1,[0 1]);title(sprintf('BM3D1: PSNR = %g',myPSNR(IO,recB1,1)));
t5 = nexttile;imagesc(recG2,[0 1]);title(sprintf('GBM3D2: PSNR = %g',myPSNR(IO,recG2,1)));
t6 = nexttile;imagesc(recB2,[0 1]);title(sprintf('BM3D2: PSNR = %g',myPSNR(IO,recB2,1)));
t7 = nexttile;imagesc(recGB,[0 1]);title(sprintf('G1-B2: PSNR = %g',myPSNR(IO,recGB,1)));
t8 = nexttile;imagesc(recBG,[0 1]);title(sprintf('B1-G2: PSNR = %g',myPSNR(IO,recBG,1)));
t1 = nexttile;imagesc(IO,[0 1]);title('original');
linkaxes([t1 t2 t3 t4 t5 t6 t7 t8]);


