clear;
d = 128; % image dim.
sZ = 16; % size of blocks
SNR = 4; % SNR of noisy image

opts.numMax =16; % max number of blocks to match for each ref. block
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
path = 'C:\Users\toby.sanders\Dropbox\archives\data\testImages\';
rng(2021);

% get image and add noise
% I = imresize(im2double(rgb2gray(imread('peppers.png'))),[d,d]);
% I = im2double(imread('cameraman.tif'));
% I = im2double(imread('images/barbara.png'));
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
gpuID = gpuDevice(2);
reset(gpuID);

tic;
[recBM4,out] = LTBM3D_distributed(I,sigma,opts);
toc;

tic;
[recBM5,out] = LTBM3D(I,sigma,opts);
toc;

tic;
recBM = BM3D(I,sigma);
toc;

myrel(recBM4,IO)
myrel(recBM5,IO)
psnr1 = myPSNR(IO,recBM4);
psnr2 = myPSNR(IO,recBM5);
psnr3 = myPSNR(IO,recBM);
%%
figure(201);colormap(gray);tiledlayout(2,2,'tilespacing','none');
t1 = nexttile;imagesc(recBM,[0 1]);title(sprintf('original code: PSNR = %g',psnr3));
t2 = nexttile;imagesc(I,[0 1]);title('noisy');
t3 = nexttile;imagesc(recBM4,[0 1]);title(sprintf('parallel ver: PSNR = %g',psnr1));
t4 = nexttile;imagesc(recBM5,[0 1]);title(sprintf('normal ver: PSNR = %g',psnr2));
linkaxes([t1 t2 t3 t4]);