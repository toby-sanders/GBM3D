clear;
d = 256; % image dim.
sZ = 32; % size of blocks
SNR = 5; % SNR of noisy image
opts.blockSize = sZ;
opts.numMax = 64; % max number of blocks to match for each ref. block
opts.numMin = 32; % min number of blocks to match for each ref. block
opts.wname = 'sym';
opts.wnamez = 'db';
opts.order = 4;
opts.orderz = 1;
opts.levels = 3;

rng(2021);

% get image and add noise
I = imresize(im2double(rgb2gray(imread('peppers.png'))),[d,d]);
% I = im2double(imread('cameraman.tif'));
% I = im2double(imread('images/barbara.png'));
% I = im2double(rgb2gray(imread('images/IMG_0086.JPG')));
% I = I(1201:1201+d-1,1601:1601+d-1);
% I = phantom(d);
% I = im2double(imread('C:\Users\toby.sanders\Desktop\cameraman.bmp'));
% I = im2double(rgb2gray(imread('images/lena.png')));
% I = im2double(rgb2gray(imread('images/peppers2.png')));
% I = im2double(rgb2gray(imread('images/monarch.png')));
% I = im2double(rgb2gray(imread('images/tulips.png')));
% I = imresize(im2double(rgb2gray(imread('images/cat.png'))),[768 512]);

IO = I;
[I,sigma] = add_Wnoise(I,SNR);

d = size(I,1);
% match blocks then denoise
tic;
[S] = matchBlocks(I,sigma,opts);



[U,Volume,out] = denoiseAndAggregate(I,S,sZ,sigma,opts);
toc;

% tic;
% S2 = matchBlocksGPU(U1,sigma/4,opts);
% toc;
% 
% [U,Volume,out] = denoiseAndAggregateGPU(I,S2,sZ,sigma,opts);

dd = squeeze(sum(sum(abs(diff(Volume,1,3)))))./squeeze(sum(sum(abs(Volume(:,:,1:end-1)))));

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

figure(111);colormap(gray);
tiledlayout(2,2,'tilespacing','none');
t1 = nexttile;imagesc(I,[0 1]);colorbar;title('noisy');
set(gca,'fontsize',20);
t2 = nexttile;imagesc(U,[0 1]);colorbar;title('my Block matching denoised');
set(gca,'fontsize',20);% imcontrast;
% t3 = nexttile;imagesc(recTV,[0 1]);colorbar;title('TV denoise');
t3 = nexttile;imagesc(IO,[0 1]);title('original');
set(gca,'fontsize',20);
t4 = nexttile;imagesc(recBM,[0 1]);colorbar;title('the real BM3D');
set(gca,'fontsize',20);
linkaxes([t1 t2 t3 t4]);


figure(112);colormap(jet);tiledlayout(2,2,'tilespacing','none');
nexttile;imagesc(out.Z);title('Z, num matches at each ref.');colorbar;
nexttile;imagesc(out.K);title('K, num matches at each pixel');colorbar;
% nexttile;imagesc(out1.Z);title('Z, num matches at each ref. (first pass)');
colorbar;
% nexttile;imagesc(out1.K);title('K, num matches at each pixel (first pass)');
nexttile;plot(out.dd);hold on;
plot(dd);hold off;
colorbar;