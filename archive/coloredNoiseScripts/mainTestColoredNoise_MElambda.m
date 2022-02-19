clear;
d = 256; % image dim.
sZ = 32; % size of blocks
SNR = 100; % SNR of noisy image
opts.blockSize = sZ;
opts.numMax = 16; % max number of blocks to match for each ref. block
opts.numMin = 16; % min number of blocks to match for each ref. block
opts.wname = 'bior';
opts.wnamez = 'db';
opts.order =13;
opts.orderz = 1;
opts.levels = 3;
opts.cycleSpin = 2;
opts.matchSpins =  [0 4 sZ/2-1];
opts.wiener = false;
opts.filtType = 'ht';
omega = 3;
lambda = 50;
rng(2021);

path = 'C:\Users\toby.sanders\Dropbox\TobySharedMATLAB\TobyShared\solvers\DenoisingEngines\LTBM3D\images\';
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

parm = dopts;
parm.theta = sigma^2*lambda;
% [U0,out0] = HOTVL2_deblur(h,b,parm);
[U0,out0] = LTBM3D_deconvME(b,hhat,sigma,dopts);
lambda0 = out0.ME.thetas(end)/out0.ME.sigmas(end)

dopts.lambda = lambda0*5;
dopts.order = 1;
[U1,out] = LTBM3D_deconv(b,hhat,sigma,dopts);
% [U1] = LTBM3D(recWei,sigmaPSD,opts);
% U2 = LTBM3D(recWei,8*mean(sigma*col(abs(filt))),opts);

U3 = BM3DDEB(b,sigma,fftshift(h));

%%
fprintf('ME soln = %g\n',myrel(U0,I));
fprintf('my 2 step soln = %g\n',myrel(U1,I));
fprintf('BM3DDEB soln = %g\n',myrel(U3,I));
fprintf('first of 2 step soln = %g\n',myrel(out.U0,I));
figure(567);colormap(gray);tiledlayout(3,3,'tilespacing','none');
t0 = nexttile;imagesc(b,[0 1]);title('blurry data');
t1 = nexttile;imagesc(out.recWei1,[0 1]);title('weiner filter soln');
t2 = nexttile;imagesc(U1,[0 1]);title('my colored BM3d');
% t3 = nexttile;imagesc(U2,[0 1]);title('my white BM3D');
t6 = nexttile;imagesc(I,[0 1]);title('original');
t5 = nexttile;imagesc(fftshift(out.sigmaPSD));colorbar;title('first PSD');
t4 = nexttile;imagesc(U3,[0 1]);title('original colored BM3DDEB');
t7 = nexttile;imagesc(out.U0,[0 1]);title('initial soln');
t8 = nexttile;imagesc(out.recWei2,[0 1]);title('second Weiner soln');
% t9 = nexttile;imagesc(out.sigmaPSD2);title('second PSD');colorbar;
t9 = nexttile;imagesc(U0,[0 1]);title('Me soln');
linkaxes([t0 t1 t2 t4 t6 t7 t8 t9]);
