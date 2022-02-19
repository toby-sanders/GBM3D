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
opts.tauMode = 2;
rng(2021);

% get image and add noise
% I = imresize(im2double(rgb2gray(imread('peppers.png'))),[d,d]);
% I = im2double(imread('cameraman.tif'));
% I = im2double(imread('images/barbara.png'));
% I = im2double(rgb2gray(imread('images/IMG_0086.JPG')));
% I = I(1201:1201+d-1,1601:1601+d-1);
% I = phantom(d);
% I = im2double(imread('C:\Users\toby.sanders\Desktop\cameraman.bmp'));
% I = im2double(rgb2gray(imread('images/lena.png')));
I = im2double(rgb2gray(imread('images/peppers2.png')));
% I = im2double(rgb2gray(imread('images/monarch.png')));
% I = im2double(rgb2gray(imread('images/tulips.png')));
% I = imresize(im2double(rgb2gray(imread('images/cat.png'))),[768 512]);
gpuID = gpuDevice(2);
reset(gpuID);

IO = I;
[I,sigma] = add_Wnoise(I,SNR);
tic;
[U1] = LTBM3D(I,sigma,opts);
toc;
tic;
[U2] = LTBM3D2(I,sigma,opts);
toc;

%%
myrel(U1,IO,1)
myrel(U2,IO,1)
figure(1212);colormap(gray);
tiledlayout(2,2,'tilespacing','none');
t1 = nexttile;imagesc(IO,[0 1]);title('original');
t2 = nexttile;imagesc(I,[0 1]);
t3 = nexttile;imagesc(U1,[0 1]);title('single filtered');
t4 = nexttile;imagesc(U2,[0 1]);title('double filtered');
linkaxes([t1 t2 t3 t4]);

