clear;
d = 128; % image dim.
SNR = 12; % SNR of noisy image
path = 'C:\Users\toby.sanders\Dropbox\archives\data\testImages\';
rng(2021);

% get image and add noise
% I = imresize(im2double(rgb2gray(imread('peppers.png'))),[d,d]);
% I = im2double(imread('cameraman.tif'));
% I = im2double(imread([path,'barbara.png']));
% I = im2double(rgb2gray(imread([path,'IMG_0086.JPG'])));
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


opts.profile = 'superFast';
tic;U1 = GBM3D_distributed(I,sigma,opts);toc;

opts.profile = 'fast';
tic;U2 = GBM3D_distributed(I,sigma,opts);toc;

opts.profile = 'default';
tic;U3 = GBM3D_distributed(I,sigma,opts);toc;

opts.profile = 'accuracy';
tic;U4 = GBM3D_distributed(I,sigma,opts);toc;


%%
figure(42);tiledlayout(2,2,'tilespacing','compact');colormap(gray);
t1 = nexttile;imagesc(U1,[0 1]);
title(sprintf('super fast, PSNR = %g',myPSNR(IO,U1,1)));
t2 = nexttile;imagesc(U2,[0 1]);
title(sprintf('fast, PSNR = %g',myPSNR(IO,U2,1)));
t3 = nexttile;imagesc(U3,[0 1]);
title(sprintf('default, PSNR = %g',myPSNR(IO,U3,1)));
t4 = nexttile;imagesc(U4,[0 1]);
title(sprintf('accuracy, PSNR = %g',myPSNR(IO,U4,1)));
linkaxes([t1 t2 t3 t4]);