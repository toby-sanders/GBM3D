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
opts.tauMode = 1;
path = 'C:\Users\toby.sanders\Dropbox\archives\data\testImages\';
rng(2021);

% get image and add noise
% I = imresize(im2double(rgb2gray(imread('peppers.png'))),[d,d]);
% I = im2double(imread('cameraman.tif'));
% I = im2double(imread('images/barbara.png'));
% I = im2double(rgb2gray(imread('images/IMG_0086.JPG')));
% I = I(1201:1201+d-1,1601:1601+d-1);
% I = phantom(d);
I = im2double(rgb2gray(imread([path,'lena.png'])));
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
[U] = GBM3D(I,sigma,opts);
toc;

% denoise for baseline comparisons
tic;
profile = 'np';
recBM = BM3D(I,sigma,profile,BM3DProfile.HARD_THRESHOLDING);
toc;
tic;
recBM2 = BM3D(I,sigma);
toc;


h = zeros(size(I,1),size(I,2));h(1) = 1;
hopts.order = 2;
hopts.levels = 1;
hopts.iter = 75;
hopts.wrap_shrink = false;
hopts.order = 1;
hopts.mode = 'deconv';
lambda = .1567; % factor for TV parameter
hopts.mu = lambda/sigma^2;
recTV = HOTV3D(h,I,[size(I,1),size(I,2),1],hopts);

%% display
fprintf('my error: %g\n',myrel(U,IO,1));
fprintf('BM3D error: %g\n',myrel(recBM,IO,1));
fprintf('BM3D 2 step %g\n',myrel(recBM2,IO,1));
fprintf('TV error: %g\n',myrel(recTV,IO,1));
% fprintf('first pass: %g\n',myrel(U1,IO,1));


psnr1 = myPSNR(IO,I,1);
psnr2 = myPSNR(IO,U,1);
psnr3 = myPSNR(IO,recBM,1);
psnr4 = myPSNR(IO,recBM2,1);
% figure(114);colormap(gray);
% tiledlayout(1,3,'tilespacing','none');
% f1 = nexttile;imagesc(I,[0 1]);title({'noisy image';sprintf('PSNR = %g',psnr1)});
% set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);
% f2 = nexttile;imagesc(U,[0 1]);title({'denoised with our prototype';sprintf('PSNR = %g',psnr2)});
% set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);
% f3 = nexttile;imagesc(recBM,[0 1]);title({'denoised with original BM3D';sprintf('PSNR = %g',psnr3)});
% set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);
% linkaxes([f1 f2 f3]);


figure(116);colormap(gray);
tiledlayout(2,2,'tilespacing','none');
f2 = nexttile;imagesc(U,[0 1]);title({'denoised with our prototype';sprintf('PSNR = %g',psnr2)});
set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);
f3 = nexttile;imagesc(recBM,[0 1]);title({'denoised with original BM3D';sprintf('PSNR = %g',psnr3)});
set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);
f1 = nexttile;imagesc(I,[0 1]);title({'noisy image';sprintf('PSNR = %g',psnr1)});
set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);
f4 = nexttile;imagesc(recBM2,[0 1]);title({'denoised with 2 step BM3D';sprintf('PSNR = %g',psnr4)});
set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);
linkaxes([f1 f2 f3 f4]);