clear;
d = 256; % image dim.
sZ = 32; % size of blocks
SNR = 4; % SNR of noisy image
opts.blockSize = sZ;
opts.numMax = 32; % max number of blocks to match for each ref. block
opts.numMin = 20; % min number of blocks to match for each ref. block
opts.wname = 'bior';
opts.wnamez = 'db';
opts.order = 13;
opts.orderz = 1;
opts.levels = 3;
opts.cycleSpin = 2;
opts.matchSpins =  [0 sZ/2-1];
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
% I = im2double(rgb2gray(imread('images/lena.png')));
% I = im2double(rgb2gray(imread('images/peppers2.png')));
% I = im2double(rgb2gray(imread('images/monarch.png')));
% I = im2double(rgb2gray(imread('images/tulips.png')));
I = im2double(rgb2gray(imread('images/cat.png')));
% I = im2double(rgb2gray(imread('images/baboon.png')));

IO = I;
[I,sigma] = add_Wnoise(I,SNR);

% our BM3D code
gpuID = gpuDevice(2);
reset(gpuID);
tic;
[U] = LTBM3D(I,sigma,opts);
toc;

% denoise for baseline comparisons
tic;
profile = 'np';
recBM = BM3D(I,sigma,profile,BM3DProfile.HARD_THRESHOLDING);
toc;

opts.wname = 'sym';
opts.order = 4;
tic;
[U2] = LTBM3D(I,sigma,opts);
toc;


%% display
fprintf('my error bior: %g\n',myrel(U,IO,1));
fprintf('my error sym %g\n',myrel(U2,IO,1));
fprintf('BM3D error: %g\n',myrel(recBM,IO,1));



psnr1 = myPSNR(IO,I,1);
psnr2 = myPSNR(IO,U,1);
psnr3 = myPSNR(IO,recBM,1);
psnr4 = myPSNR(IO,U2,1);


figure(116);colormap(gray);
tiledlayout(2,2,'tilespacing','none');
f2 = nexttile;imagesc(U,[0 1]);title({'biorthogonal';sprintf('PSNR = %g',psnr2)});
set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);
f4 = nexttile;imagesc(U2,[0 1]);title({'sym wavelets';sprintf('PSNR = %g',psnr4)});
set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);
f3 = nexttile;imagesc(recBM,[0 1]);title({'denoised with original BM3D';sprintf('PSNR = %g',psnr3)});
set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);
f1 = nexttile;imagesc(I,[0 1]);title({'noisy image';sprintf('PSNR = %g',psnr1)});
set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);

linkaxes([f1 f2 f3 f4]);