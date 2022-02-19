clear;
d = 256; % image dim.
sZ = 32; % size of blocks
SNR = 4; % SNR of noisy image
opts.blockSize = sZ;
opts.numMax = 16; % max number of blocks to match for each ref. block
opts.numMin = 16; % min number of blocks to match for each ref. block
opts.wname = 'bior';
opts.wnamez = 'db';
opts.order = 13;
opts.orderz = 1;
opts.levels = 3;
opts.cycleSpin = 2;
opts.matchSpins =  [0 sZ/2-1];
opts.wiener = false;
rng(2021);

% get image and add noise
% I = imresize(im2double(rgb2gray(imread('peppers.png'))),[d,d]);
% I = im2double(imread('cameraman.tif'));
% I = im2double(imread('images/barbara.png'));
I = im2double(rgb2gray(imread('images/IMG_0086.JPG')));
I = I(1201:1201+d-1,1601:1601+d-1);
% I = phantom(d);
% I = im2double(rgb2gray(imread('images/lena.png')));
% I = im2double(rgb2gray(imread('images/peppers2.png')));
% I = im2double(rgb2gray(imread('images/monarch.png')));
% I = im2double(rgb2gray(imread('images/tulips.png')));
% I = im2double(rgb2gray(imread('images/cat.png')));
% I = im2double(rgb2gray(imread('images/baboon.png')));

IO = I;


% our BM3D code
gpuID = gpuDevice(2);
reset(gpuID);

ee1 = zeros(4,1);
SNR = 4;
[I,sigma] = add_Wnoise(IO,SNR);
opts.tauMode = 1;
[U1] = LTBM3D(I,sigma,opts);
ee1(1) = myrel(U1,IO);
opts.tauMode = 2;
[U2] = LTBM3D(I,sigma,opts);
ee1(2) = myrel(U2,IO);
opts.tauMode = 3;
[U3] = LTBM3D(I,sigma,opts);
ee1(3) = myrel(U3,IO);
opts.tauMode = 4;
[U4] = LTBM3D(I,sigma,opts);
ee1(4) = myrel(U4,IO);


ee2 = ee1;
SNR = 10;
[I,sigma] = add_Wnoise(IO,SNR);
opts.tauMode = 1;
[U5] = LTBM3D(I,sigma,opts);
ee2(1) = myrel(U5,IO);
opts.tauMode = 2;
[U6] = LTBM3D(I,sigma,opts);
ee2(2) = myrel(U6,IO);
opts.tauMode = 3;
[U7] = LTBM3D(I,sigma,opts);
ee2(3) = myrel(U7,IO);
opts.tauMode = 4;
[U8] = LTBM3D(I,sigma,opts);
ee2(4) = myrel(U8,IO);


%%
figure(123);colormap(gray);
tiledlayout(2,4,'tilespacing','none');
t1 = nexttile;imagesc(U1,[0 1]);
t2 = nexttile;imagesc(U2,[0 1]);
t3 = nexttile;imagesc(U3,[0 1]);
t4 = nexttile;imagesc(U4,[0 1]);
t5 = nexttile;imagesc(U5,[0 1]);
t6 = nexttile;imagesc(U6,[0 1]);
t7 = nexttile;imagesc(U7,[0 1]);
t8 = nexttile;imagesc(U8,[0 1]);
linkaxes([t1 t2 t3 t4 t5 t6 t7 t8]);

figure(124);
subplot(1,2,1);plot(ee1);
subplot(1,2,2);plot(ee2);
ee1
ee2