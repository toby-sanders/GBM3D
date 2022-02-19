clear;
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
opts.cycleSpin = 2;
opts.matchSpins = [0 sZ/2-1];
opts.wiener = false;
rng(2021);

% get image

I = im2double(imread('images/soccerField.tif'));
% I = I(:,1:512);
% I = imresize(I,[512,512]);
sigma = .025;

I = im2double(imread('images/soccerField2.tif'));
% I = I(:,1:512);
% I = imresize(I,[512,512]);
sigma = .017;


% I = im2double(imread('images/manOnRoof.tif'));
% I = I(1:512,1:1024);
% sigma = .055;
% 
% I = im2double(imread('images/manOnRoof2.tif'));
% I = I(1:512,1:1024);
% sigma = .045;
% 

% match blocks then denoise

tic;
[U,Volume,out] = LTBM3D(I,sigma,opts);
toc;


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
mm1 = .40;mm2 = .75;
% mm1 = 0;mm2 = 1;
figure(111);colormap(gray);
tiledlayout(2,2,'tilespacing','none');
t1 = nexttile;imagesc(I,[mm1 mm2]);colorbar;title('noisy');
set(gca,'fontsize',20);
t2 = nexttile;imagesc(U,[mm1 mm2]);colorbar;title('my Block matching denoised');
set(gca,'fontsize',20);% imcontrast;
% t3 = nexttile;imagesc(recTV,[0 1]);colorbar;title('TV denoise');
% t3 = nexttile;imagesc(IO,[0 1]);title('original');
% set(gca,'fontsize',20);
t4 = nexttile;imagesc(recBM,[mm1 mm2]);colorbar;title('the real BM3D');
set(gca,'fontsize',20);
t3 = nexttile;imagesc(recBM2,[mm1 mm2]);title('2 step BM3D');
set(gca,'fontsize',20);
linkaxes([t1 t2 t3 t4]);


figure(114);colormap(gray);
tiledlayout(2,3,'tilespacing','none');
s1 = nexttile;imagesc(I,[mm1,mm2]);title('original noisy frame');
set(gca,'fontsize',20,'Xtick',[],'Ytick',[]);
s2 = nexttile;imagesc(recBM,[mm1,mm2]);title('denoised with black box code');
set(gca,'fontsize',20,'Xtick',[],'Ytick',[]);
s3 = nexttile;imagesc(U,[mm1,mm2]);title('denoised with our prototype');
set(gca,'fontsize',20,'Xtick',[],'Ytick',[]);
linkaxes([s1,s2,s3]);

s4 = nexttile;imagesc(I,[mm1,mm2]);title('zoom of above');
set(gca,'fontsize',20,'Xtick',[],'Ytick',[]);
s5 = nexttile;imagesc(recBM,[mm1,mm2]);title('zoom of above');
set(gca,'fontsize',20,'Xtick',[],'Ytick',[]);
s6 = nexttile;imagesc(U,[mm1,mm2]);title('zoom of above');
set(gca,'fontsize',20,'Xtick',[],'Ytick',[]);
linkaxes([s4 s5 s6]);