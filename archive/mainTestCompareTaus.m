clear;
d = 256; % image dim.
sZ = 32; % size of blocks
SNR = 10; % SNR of noisy image
opts.blockSize = sZ;
opts.numMax = 32; % max number of blocks to match for each ref. block
opts.numMin = 20; % min number of blocks to match for each ref. block
opts.wname = 'sym';
opts.wnamez = 'db';
opts.order = 4;
opts.orderz = 1;
opts.levels = 3;
opts.cycleSpin = 2;
opts.matchSpins =  [0 sZ/2-1];
opts.wiener = false;
rng(2021);

% get image and add noise
% I = im2double(rgb2gray(imread('peppers.png')));
% I = im2double(imread('cameraman.tif'));
% I = im2double(imread('images/barbara.png'));
% I = im2double(rgb2gray(imread('images/IMG_0086.JPG')));
% I = I(1201:1201+d-1,1601:1601+d-1);
% I = phantom(d);
% I = im2double(rgb2gray(imread('images/lena.png')));
% I = im2double(rgb2gray(imread('images/peppers2.png')));
I = im2double(rgb2gray(imread('images/monarch.png')));
% I = im2double(rgb2gray(imread('images/tulips.png')));
% I = im2double(rgb2gray(imread('images/cat.png')));
% I = im2double(rgb2gray(imread('images/baboon.png')));

IO = I;
[I,sigma] = add_Wnoise(I,SNR);

% our BM3D code
gpuID = gpuDevice(2);
reset(gpuID);

opts.tauMode = 1;
[U1] = LTBM3D(I,sigma,opts);
opts.tauMode = 2;
[U2] = LTBM3D(I,sigma,opts);
opts.tauMode = 3;
[U3] = LTBM3D(I,sigma,opts);

%% display
fprintf('case 1 error: %g\n',myrel(U1,IO,1));
fprintf('case 2 error: %g\n',myrel(U2,IO,1));
fprintf('case 3 error %g\n',myrel(U3,IO,1));



psnr1 = myPSNR(IO,I,1);
psnr2 = myPSNR(IO,U1,1);
psnr3 = myPSNR(IO,U2,1);
psnr4 = myPSNR(IO,U3,1);



figure(116);colormap(gray);
tiledlayout(2,2,'tilespacing','none');
f2 = nexttile;imagesc(U1,[0 1]);title({'case 1';sprintf('PSNR = %g',psnr2)});
set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);
f3 = nexttile;imagesc(U2,[0 1]);title({'case 2';sprintf('PSNR = %g',psnr3)});
set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);
f1 = nexttile;imagesc(I,[0 1]);title({'noisy image';sprintf('PSNR = %g',psnr1)});
set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);
f4 = nexttile;imagesc(U3,[0 1]);title({'case 3';sprintf('PSNR = %g',psnr4)});
set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);
linkaxes([f1 f2 f3 f4]);