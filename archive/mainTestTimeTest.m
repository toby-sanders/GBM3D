clear;
d = 512; % image dim.
sZ = 32; % size of blocks
SNR = 4; % SNR of noisy image
opts.blockSize = sZ;
opts.numMax = 16; % max number of blocks to match for each ref. block
opts.numMin = 16; % min number of blocks to match for each ref. block
opts.wname = 'sym';
opts.wnamez = 'db';
opts.order = 4;
opts.orderz = 1;
opts.levels = 3;
opts.cycleSpin = 2;
opts.matchSpins = [0];
opts.wiener = false;
rng(2021);

% get image and add noise
% I = imresize(im2double(rgb2gray(imread('peppers.png'))),[d,d]);
% I = im2double(imread('cameraman.tif'));
% I = im2double(imread('images/barbara.png'));
% I = im2double(rgb2gray(imread('images/IMG_0086.JPG')));
% I = I(1201:1201+d-1,1601:1601+d-1);
% I = phantom(d);
% I = im2double(imread('C:\Users\toby.sanders\Desktop\cameraman.bmp'));
II = im2double(rgb2gray(imread('images/lena.png')));
% I = im2double(rgb2gray(imread('images/peppers2.png')));
% I = im2double(rgb2gray(imread('images/monarch.png')));
% I = im2double(rgb2gray(imread('images/tulips.png')));
% I = imresize(im2double(rgb2gray(imread('images/cat.png'))),[768 512]);
dAll = [128,192, 256, 384 512 768 1024];
tt = zeros(numel(dAll),3);

for i = 1:numel(dAll)
    
IO = imresize(II,[dAll(i),dAll(i)]);
[I,sigma] = add_Wnoise(IO,SNR);
shift = sZ/2-1;
d = size(I,1);
% shifts = [0 sZ/2-1];

% our BM3D code
gpuID = gpuDevice(2);
reset(gpuID);


tic;
U = LTBM3D(I,sigma,opts);
tt(i,1) = toc;

% denoise for baseline comparisons
tic;
profile = 'np';
recBM = BM3D(I,sigma,profile,BM3DProfile.HARD_THRESHOLDING);
tt(i,2) = toc;

tic;
recBM2 = BM3D(I,sigma);
tt(i,3) = toc;
end

%%
figure(11);
subplot(1,2,1);hold off;
plot(log(dAll)/log(2),tt(:,1),'--xk','linewidth',2);hold on;
plot(log(dAll)/log(2),tt(:,2),':k*','linewidth',2);
plot(log(dAll)/log(2),tt(:,3),'-ok','linewidth',2);
xlabel('N (image dimension 2^N)');

set(gca,'Xtick',[7 8 9 10],'fontweight','bold','fontsize',16,'Ytick',[10 20 30 40]);
ylabel('execution time (s)');
legend({'our GPU variation','1 step BM3D','2 step BM3D'},'location','northwest');
title('comparison of exectution time');

subplot(1,2,2);hold off;
plot(log(dAll)/log(2),tt(:,2)./tt(:,1),'--xk','linewidth',2);hold on;
plot(log(dAll)/log(2),tt(:,3)./tt(:,1),':k*','linewidth',2);
xlabel('N (image dimension 2^N)');
axis([7 10 0 12])
set(gca,'Xtick',[7 8 9 10],'fontweight','bold','fontsize',16);
ylabel('speed up factor');
legend({'1 step BM3D','2 step BM3D'},'location','northwest');
title('speed up factor of our GPU variation');





