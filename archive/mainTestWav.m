% prototyping BM3D using wavelets
clear;
d = 512; % image dim.
sZ = 8; % size of blocks
SNR = 5; % SNR of noisy image
opts.blockSize = sZ;
opts.numMax = 256; % max number of blocks to match for each ref. block
opts.numMin = 20; % min number of blocks to match for each ref. block
my2Dtau = 0.5; % shrinkage coef. for init 2D denoise
rng(2021);

% get image and add noise
% I = imresize(im2double(rgb2gray(imread('peppers.png'))),[d,d]);
% I = im2double(imread('cameraman.tif'));
I = im2double(imread('images/barbara.png'));
% I = im2double(rgb2gray(imread('IMG_0086.JPG')));
% I = I(1201:1201+d-1,1601:1601+d-1);
% I = phantom(d);
IO = I;
[I,sigma] = add_Wnoise(I,SNR);

% initial denoise before block matching
[dec,app] = wavedec2(I,4,'sym4');
dec = dec.*(abs(dec)>my2Dtau*sigma);
Id0 = waverec2(dec,app,'sym4');

% first block matching and denoising
S = matchBlocks(Id0,sigma,opts);
dopts.numMax = opts.numMax;
[U0,out] = denoiseBlocksWav(I,S,sZ,sigma,dopts);
U0(isnan(U0)) = 1;

% re-evaluate block matching and denoise
S = matchBlocks(U0,sigma/2,opts);
dopts.numMax = opts.numMax;
[U,out] = denoiseBlocksWav(I,S,sZ,sigma,dopts);


% denoise for baseline comparisons
recBM = BM3D(I,sigma);
h = zeros(d,d);h(1) = 1;
hopts.order = 2;
hopts.levels = 1;
hopts.iter = 75;
hopts.wrap_shrink = false;
hopts.order = 1;
hopts.mode = 'deconv';
lambda = .0367; % factor for TV parameter
hopts.mu = lambda/sigma^2;
recTV = HOTV3D(h,I,[d,d,1],hopts);


%% display
zz1 = U(50:200,50:200);
zz2 = IO(50:200,50:200);
zz3 = U0(50:200,50:200);
zz4 = Id0(50:200,50:200);
zz5 = recBM(50:200,50:200);
fprintf('my second error: %g\n',norm(zz1(:)-zz2(:),1)/norm(zz2(:),1));
fprintf('my first error: %g\n',norm(zz3(:)-zz2(:),1)/norm(zz2(:),1));
fprintf('real BM3D error: %g\n',myrel(zz5,zz2,1));
% fprintf('BM3D error: %g\n',myrel(recBM(50:200,50:200),IO(50:200,50:200),1));
% fprintf('TV error: %g\n',myrel(recTV(50:200,50:200),IO(50:200,50:200),1));


figure(111);colormap(gray);
tiledlayout(2,2,'tilespacing','none');
t1 = nexttile;imagesc(I,[0 1]);colorbar;title('noisy');
set(gca,'fontsize',20);
t2 = nexttile;imagesc(U,[0 1]);colorbar;title('my Block matching denoised');
set(gca,'fontsize',20);% imcontrast;
t3 = nexttile;imagesc(U0,[0 1]);colorbar;title('first estimate');
set(gca,'fontsize',20);
t4 = nexttile;imagesc(recBM,[0 1]);colorbar;title('the real BM3D');
set(gca,'fontsize',20);
linkaxes([t1 t2 t3 t4]);


figure(112);colormap(jet);tiledlayout(2,2,'tilespacing','none');
nexttile;imagesc(out.Z);title('Z, num matches at each ref.');colorbar;
nexttile;imagesc(out.K);title('K, num matches at each pixel');colorbar;
