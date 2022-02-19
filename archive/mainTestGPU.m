clear;
d = 256; % image dim.
sZ = 8; % size of blocks
SNR = 5; % SNR of noisy image
k = 2;
opts.blockSize = sZ;
opts.numMax = 128; % max number of blocks to match for each ref. block
opts.numMin = 16; % min number of blocks to match for each ref. block
lambda = .1567; % factor for TV parameter
rng(2021);

% get image and add noise
% I = imresize(im2double(rgb2gray(imread('peppers.png'))),[d,d]);
I = im2double(imread('cameraman.tif'));
% I = im2double(imread('images/barbara.png'));
% I = im2double(rgb2gray(imread('IMG_0086.JPG')));
% I = I(1201:1201+d-1,1601:1601+d-1);
% I = phantom(d);
% I = im2double(imread('C:\Users\toby.sanders\Desktop\cameraman.bmp'));
IO = I;
[I,sigma] = add_Wnoise(I,SNR);

% match blocks then denoise
% tic;
% S = matchBlocks(I,sigma,opts);
% toc;
% dopts.numMax = opts.numMax;
% [U,out] = denoiseBlocks_GPU(I,S,sZ,sigma,dopts);
[U,out] = LTBM3D(I,sigma,opts);

% denoise for baseline comparisons
tic;
recBM = BM3D(I,sigma);
toc;
h = zeros(d,d);h(1) = 1;
hopts.order = 2;
hopts.levels = 1;
hopts.iter = 75;
hopts.wrap_shrink = false;
hopts.order = 1;
hopts.mode = 'deconv';
hopts.mu = lambda/sigma^2;
recTV = HOTV3D(h,I,[d,d,1],hopts);

%% display
fprintf('my error: %g\n',myrel(U(50:200,50:200),IO(50:200,50:200),1));
fprintf('BM3D error: %g\n',myrel(recBM(50:200,50:200),IO(50:200,50:200),1));
fprintf('TV error: %g\n',myrel(recTV(50:200,50:200),IO(50:200,50:200),1));


figure(111);colormap(gray);
tiledlayout(2,2,'tilespacing','none');
t1 = nexttile;imagesc(I,[0 1]);colorbar;title('noisy');
set(gca,'fontsize',20);
t2 = nexttile;imagesc(U,[0 1]);colorbar;title('my Block matching denoised');
set(gca,'fontsize',20);% imcontrast;
t3 = nexttile;imagesc(recTV,[0 1]);colorbar;title('TV denoise');
set(gca,'fontsize',20);
t4 = nexttile;imagesc(recBM,[0 1]);colorbar;title('the real BM3D');
set(gca,'fontsize',20);
linkaxes([t1 t2 t3 t4]);


figure(112);colormap(jet);tiledlayout(2,2,'tilespacing','none');
nexttile;imagesc(out.Z);title('Z, num matches at each ref.');
nexttile;imagesc(out.K);title('K, num matches at each pixel');