clear;
d = 128; % image dim.
sZ = 32; % size of blocks
SNR = 5; % SNR of noisy image
opts.blockSize = sZ;
opts.numMax =16; % max number of blocks to match for each ref. block
opts.numMin = 16; % min number of blocks to match for each ref. block
opts.wname = 'bior';
opts.wnamez = 'db';
opts.order = 15;
opts.orderz = 1;
opts.levels = 3;
opts.cycleSpin = 2;
opts.matchSpins =  [0 4 6];
opts.tauMode = 2;
opts.blockSize = 16;
opts.blockSizeWie = 8;
opts.matchSize = 32;
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
[U] = LTBM3D(I,sigma,opts);
toc;


% I put a break point in before the second Wiener filter to so that I can
% go in and cheat with for the first estimate by loading in the original
% image
U2 = LTBM3D(I,sigma,opts);



%%
myrel(U,IO)
myrel(U2,IO)


figure(118);tiledlayout(2,2,'tilespacing','none');colormap(gray);
t1 = nexttile;imagesc(U,[0 1]);title('2 step BM3D');
t2 = nexttile;imagesc(U2,[0 1]);title('cheated Wiener filter');
t3 = nexttile;imagesc(I,[0 1]);title('noisy');
linkaxes([t1 t2 t3]);
