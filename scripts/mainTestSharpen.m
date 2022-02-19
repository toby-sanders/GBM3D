clear;
d = 256; % image dim.
sZ = 32; % size of blocks
SNR = 55; % SNR of noisy image
opts.profile = 'default';
omega = 1.5;
lambda = 50;
alpha = 1.5;
rng(2021);

% get image and add noise
path = 'C:\Users\toby.sanders\Dropbox\archives\data\testImages\';
% get image and add noise
% I = imresize(im2double(rgb2gray(imread('peppers.png'))),[d,d]);
% I = im2double(imread('cameraman.tif'));
% I = im2double(imread('images/barbara.png'));
% I = im2double(rgb2gray(imread('images/IMG_0086.JPG')));
% I = I(1201:1201+d-1,1601:1601+d-1);
% I = phantom(d);
% I = im2double(imread('C:\Users\toby.sanders\Desktop\cameraman.bmp'));
% I = im2double(rgb2gray(imread([path,'lena.png'])));
I = im2double((imread([path,'barbara.png'])));
% I = im2double(rgb2gray(imread([path,'cat.png'])));
% I = I(256-127:256+128,256-127:256+128);
% I = im2double(rgb2gray(imread([path,'peppers2.png'])));
% I = im2double(rgb2gray(imread('images/monarch.png')));
% I = I(1:d,141:140+d);
% I = im2double(rgb2gray(imread('images/tulips.png')));
[d1,d2] = size(I);
% [b,sigma] = add_Wnoise(I,SNR);
sigma = 2/255;

U1 = LTsharpen(I,alpha,sigma);
U2 = LTsharpen2(I,alpha,sigma);


%%

myrel(U1,U2)
figure(567);colormap(gray);tiledlayout(2,2,'tilespacing','none');
% t0 = nexttile;imagesc(b,[0 1]);title('blurry data');
t1 = nexttile;imagesc(I,[0 1]);title('original');
t2 = nexttile;imagesc(U1,[0 1]);title('sharpened 1');
t3 = nexttile;imagesc(U2,[0 1]);title('sharpened 2');
linkaxes([t1 t2 t3]);