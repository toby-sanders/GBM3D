clear;
d =256; % image dim.
sZ = 32; % size of blocks
SNR = 10; % SNR of noisy image
opts.blockSize = sZ;
opts.numMax = 16; % max number of blocks to match for each ref. block
opts.numMin = 16; % min number of blocks to match for each ref. block
opts.wname = 'bior';
opts.wnamez = 'db';
opts.order =13;
opts.orderz = 1;
opts.levels = 3;
opts.cycleSpin = 2;
opts.matchSpins =  [0 4 sZ/2];
opts.wiener = false;
opts.filtType = 'ht';
omega = .01;
lambda = 1e-5;
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
% I = im2double(rgb2gray(imread('images/peppers2.png')));
I = im2double(rgb2gray(imread('images/monarch.png')));
I = I(101:100+d,141:140+d);
% I = im2double(rgb2gray(imread('images/tulips.png')));
% I = imresize(im2double(rgb2gray(imread('images/cat.png'))),[768 512]);


[h,hhat] = makeGausPSF(d,omega);
[b,sigma] = add_Wnoise(I,SNR);

V = my_Fourier_filters(2,2,d,d,1);
filt = conj(hhat)./(conj(hhat).*hhat + sigma^2*lambda*V);
recWei = real(ifft2(filt.*fft2(b)));
sigmaPSD = sigma^2*abs(filt).^2;

[S,Volume] = matchBlocksCPU(recWei,sigma,opts);
[U1] = denoiseAndAggregateGPUColored(Volume,recWei,S,sigmaPSD,opts);


[U2] = LTBM3D(b,sigma,opts);

[U3] = LTBM3D(b,sigmaPSD,opts);
%%
myrel(U3,U2)
myrel(U2,I)
myrel(U3,I)

%%
figure(111);tiledlayout(2,2,'tilespacing','none');colormap(gray);
t1 = nexttile;imagesc(recWei,[0 1]);
t2 = nexttile;imagesc(U1,[0 1]);
t3 = nexttile;imagesc(U2,[0 1]);
t4 = nexttile;imagesc(U3,[0 1]);
linkaxes([t1 t2 t3 t4]);