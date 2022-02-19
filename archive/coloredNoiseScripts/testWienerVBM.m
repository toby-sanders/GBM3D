clear;
d = 256; % image dim.
sZ = 32; % size of blocks
SNR = 10; % SNR of noisy image
opts.blockSize = sZ;
opts.numMax = 16; % max number of blocks to match for each ref. block
opts.numMin = 16; % min number of blocks to match for each ref. block
opts.wname = 'bior';
opts.wnamez = 'db';
opts.order =13;
opts.orderz = 1;
opts.levels = 2;
opts.cycleSpin = 2;
opts.matchSpins =  [0 4 sZ/2-1];
opts.wiener = false;
opts.filtType = 'ht';
omega = 3;
rng(2021);

% lambdas = [10 30 50 100];
path = 'C:\Users\toby.sanders\Dropbox\TobySharedMATLAB\TobyShared\solvers\DenoisingEngines\LTBM3D\images\';
nL = 10;
Cs = 4;
% get image and add noise
% I = imresize(im2double(rgb2gray(imread('peppers.png'))),[d,d]);
% I = im2double(imread('cameraman.tif'));
% I = im2double(imread('images/barbara.png'));
% I = im2double(rgb2gray(imread('images/IMG_0086.JPG')));
% I = I(1201:1201+d-1,1601:1601+d-1);
% I = phantom(d);
% I = im2double(imread('C:\Users\toby.sanders\Desktop\cameraman.bmp'));
% I = im2double(rgb2gray(imread([path,'lena.png'])));
% I = I(256-127:256+128,256-127:256+128);
I = im2double(rgb2gray(imread([path,'peppers2.png'])));
% I = im2double(rgb2gray(imread('images/monarch.png')));

% I = im2double(rgb2gray(imread('images/tulips.png')));
% I = imresize(im2double(rgb2gray(imread('images/cat.png'))),[768 512]);
% I = I(1:d,141:140+d);

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
ee = zeros(nL,2);

parm = dopts;
parm.theta = sigma^2*100;
[U0,out0] = HOTVL2_deblur(h,b,parm);
lambda0 = out0.thetas(end)/out0.sigmas(end);
lambdas = linspace(lambda0/10,lambda0*1.5,nL);
mmax = 0;
mmin = 1;
for i = 1:nL
    dopts.C = Cs(1);
    dopts.lambda = lambdas(i);
    [U,out] = LTBM3D_deconv1(b,hhat,sigma,dopts);
    ee(i,1) = myrel(U,I);
    ee(i,2) = myrel(out.recWei1,I);
    if ee(i,1)>mmax
        Uworst = U;
        mmax = ee(i,1);
    end
    if ee(i,1)<mmin
        Ubest = U;
        mmin = ee(i,1);
    end
end

%%
[ppV1,pp1] = min(ee(:,2));
[ppV,pp] = min(ee(:,1));
fprintf('best L2 soln = %g\n',ppV1);
fprintf('best BM3DDEB soln = %g\n',ppV);
fprintf('corresponding L2 soln = %g\n',ee(pp,2));
fprintf('best L2 lambda = %g\n',lambdas(pp1));
fprintf('best BM3DDEB lambda = %g\n',lambdas(pp));
fprintf('ratio of lambdas = %g\n',lambdas(pp1)/lambdas(pp));
figure(111);hold off;
plot(lambdas,ee(:,1));hold on;
plot(lambdas,ee(:,2));
hold off;

figure(112);colormap(gray);
tiledlayout(2,2,'tilespacing','none');
t1 = nexttile;imagesc(Uworst,[0 1]);title('worst soln');
t2 = nexttile;imagesc(Ubest,[0 1]);title('best soln');
t3 = nexttile;imagesc(U0,[0 1]);title('ME soln');
t4 = nexttile;imagesc(b);title('data');
    
