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
opts.tauMode = 1;
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



% denoise for baseline comparisons
tic;
profile = 'np';
recBM = BM3D(I,sigma,profile,BM3DProfile.HARD_THRESHOLDING);
toc;


tic;
recBM3 = BM3D(I,sigma);
toc;


%% display
opts.Wiener = true;

[recBM4,out] = LTBM3D2(I,sigma,opts);



[recBM5,out] = LTBM3D(I,sigma,opts);



fprintf('my error: %g\n',myrel(U,IO,1));
fprintf('BM3D1 error: %g\n',myrel(recBM,IO,1));
% fprintf('our BM3D copy: %g\n',myrel(recBM2,IO,1));
fprintf('BM3D2 error: %g\n',myrel(recBM3,IO,1));
fprintf('BM3D hybrid: %g\n',myrel(recBM4,IO,1));
fprintf('BM3D hybrid2: %g\n',myrel(recBM5,IO,1));

psnr1 = myPSNR(IO,I,1);
psnr2 = myPSNR(IO,U,1);
psnr3 = myPSNR(IO,recBM,1);
psnr4 = myPSNR(IO,recBM4,1);
psnr5 = myPSNR(IO,recBM3,1);
psnr6 = myPSNR(IO,recBM5,1);


figure(116);colormap(gray);
tiledlayout(2,3,'tilespacing','none');
f2 = nexttile;imagesc(U,[0 1]);title({'denoised with our prototype';sprintf('PSNR = %g',psnr2)});
set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);
f3 = nexttile;imagesc(recBM,[0 1]);title({'denoised with original BM3D';sprintf('PSNR = %g',psnr3)});
set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);
f1 = nexttile;imagesc(I,[0 1]);title({'noisy image';sprintf('PSNR = %g',psnr1)});
set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);
f4 = nexttile;imagesc(recBM4,[0 1]);title({'parallel LTBM3D';sprintf('PSNR = %g',psnr4)});
set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);
f6 = nexttile;imagesc(recBM5,[0 1]);title({'LTBM3D2';sprintf('PSNR = %g',psnr6)});
set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);
f5 = nexttile;imagesc(recBM3,[0 1]);title({'BM3D2';sprintf('PSNR = %g',psnr5)});
set(gca,'fontsize',16,'Xtick',[],'Ytick',[]);

linkaxes([f1 f2 f3 f4 f5 f6]);