clear;
d = 256; % image dim.
sZ = 32; % size of blocks
SNR = 4; % SNR of noisy image
opts.profile = 'default';
omega = 2.5;
lambda = 50;
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
I = im2double(rgb2gray(imread([path,'lena.png'])));
% I = I(256-127:256+128,256-127:256+128);
% I = im2double(rgb2gray(imread([path,'peppers2.png'])));
% I = im2double(rgb2gray(imread('images/monarch.png')));
% I = I(1:d,141:140+d);
% I = im2double(rgb2gray(imread('images/tulips.png')));
% I = imresize(im2double(rgb2gray(imread('images/cat.png'))),[768 512]);
% I = imresize(I,[256 256]);
[d1,d2] = size(I);
[b,sigma] = add_Wnoise(I,SNR);
n = b-I;
[g,ghat] = makeGausPSF([d1,d2],1);
% xx = linspace(-d2/2,d2/2,d2);
% yy = linspace(-d1/2,d1/2,d1);
% [X,Y] = meshgrid(xx,yy);
% r1 = 100;
% R2 = 200;
% % ghat = abs(X).^2 + abs(Y).^2 > r1^2 & abs(X).^2 + abs(Y).^2 < R2^2;
% ghat = min(1./(abs(X/100)+abs(Y/100)),1);
% ghat = fftshift(double(ghat));
eta = real(ifft2(fft2(n).*ghat));
b = I + eta;

sigmaPSD = sigma^2*abs(ghat).^2;
[U1,out1] = GBM3D(b,sigmaPSD,opts);
[U2,out2] = BM3D(b,sigmaPSD*d1*d2);
[U3,out3] = BM3D(b,sqrt(mean(abs(eta(:)).^2)));

%%
myrel(I,U1)
myrel(I,U2)
myrel(I,U3)

figure(72);tiledlayout(2,3,'tilespacing','none');colormap(gray);
t1 = nexttile;imagesc(I,[0 1]);
t2 = nexttile;imagesc(b,[0 1]);
t5 = nexttile;imagesc(U3,[0 1]);
t3 = nexttile;imagesc(U1,[0 1]);
t4 = nexttile;imagesc(U2,[0 1]);
t6 = nexttile;imagesc(fftshift(ghat));
linkaxes([t1 t2 t3 t4 t5])