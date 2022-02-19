clear;
d = 256; % image dim.
sZ = 32; % size of blocks
SNR = 10; % SNR of noisy image
opts.profile = 'default';
omega = 1.5;
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
[h,hhat] = makeGausPSF([d1,d2],omega);
b = ifft2(fft2(I).*hhat);
[b,sigma] = add_Wnoise(b,SNR);

% V = my_Fourier_filters(1,1,d1,d2,1);
% filt = conj(hhat)./(conj(hhat).*hhat + sigma^2*lambda*V);
% recWei = real(ifft2(filt.*fft2(b)));
% sigmaPSD = sigma^2*abs(filt).^2;
dopts.order = 1;
dopts.levels = 1;
dopts.lambda = lambda;
dopts.C = 6;
% [U1,out] = LTBM3D_deconv(b,hhat,sigma,dopts);
[U2,out2] = GBM3D_deconv(b,hhat,sigma,opts);
U3 = BM3DDEB(b,sigma,fftshift(h));


%%
% myrel(U1,I)
myrel(U2,I)
myrel(U3,I)

% psnr1 = myPSNR(I,U1,1);
psnr2 = myPSNR(I,U2,1);
psnr3 = myPSNR(I,U3,1);

figure(567);colormap(gray);tiledlayout(3,3,'tilespacing','none');
t0 = nexttile;imagesc(b,[0 1]);title('blurry data');
t1 = nexttile;imagesc(out2.recWie,[0 1]);title('1st wiener filter soln');
t17 = nexttile;imagesc(out2.recWie2,[0 1]);title('2nd wiener solution');
% t2 = nexttile;imagesc(U1,[0 1]);title(sprintf('my 1 step BM3d: PSNR = %g',psnr1));
t9 = nexttile;imagesc(U2,[0 1]);title(sprintf('my 2 step GBM3D: PSNR = %g',psnr2));
% t5 = nexttile;imagesc(fftshift(out.sigmaPSD));colorbar;title('first PSD');
t4 = nexttile;imagesc(U3,[0 1]);title(sprintf('original BM3DDEB: PSNR = %g',psnr3));
t6 = nexttile;imagesc(I,[0 1]);title('original');
t65 = nexttile;imagesc(out2.U1,[0 1]);title(sprintf('1st soln: PSNR = %g',myPSNR(I,out2.U1,1)));
% t7 = nexttile;imagesc(out.U0,[0 1]);title('initial soln');
% t9 = nexttile;imagesc(out.sigmaPSD);title('second PSD');colorbar;
linkaxes([t0 t1  t4 t6 t17 t9]);
