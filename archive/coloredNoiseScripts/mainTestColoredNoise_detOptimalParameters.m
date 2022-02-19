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
omega = 2;
rng(2021);

lambdas = [10 30 50 100];
Cs = [3 4 5];
% get image and add noise
% I = imresize(im2double(rgb2gray(imread('peppers.png'))),[d,d]);
% I = im2double(imread('cameraman.tif'));
% I = im2double(imread('images/barbara.png'));
% I = im2double(rgb2gray(imread('images/IMG_0086.JPG')));
% I = I(1201:1201+d-1,1601:1601+d-1);
% I = phantom(d);
% I = im2double(imread('C:\Users\toby.sanders\Desktop\cameraman.bmp'));
I = im2double(rgb2gray(imread('images/lena.png')));
I = I(256-127:256+128,256-127:256+128);
% I = im2double(rgb2gray(imread('images/peppers2.png')));
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
ee = zeros(numel(Cs),numel(lambdas));
mmax = 0;
mmin = 1;
for i = 1:numel(Cs)
    for j = 1:numel(lambdas)
        dopts.lambda = lambdas(j);
        dopts.C = Cs(i);
        [U1,out] = LTBM3D_deconv(b,hhat,sigma,dopts);
        ee(i,j) = myrel(U1,I);
        if ee(i,j)>mmax
            mmax = ee(i,j);
            Uworst = U1;
        end
        if ee(i,j)<mmin
            mmin = ee(i,j);
            Ubest = U1;
        end
            
    end
    Ufinal = U1;
end

U2 = BM3DDEB(b,sigma,fftshift(h));
eeO = myrel(U2,I);

sname = sprintf('errorsMonarch_SNR=%i_omega=%i',SNR,omega);
cd ColoredResults;
save(sname,'ee','omega','SNR','eeO','lambdas','Cs');
cd ../
%%
figure(687);tiledlayout(2,2,'tilespacing','none');
colormap(gray);
t1 = nexttile;
imagesc(b,[0 1]);title('data');
t2 = nexttile;
imagesc(Uworst);title('worst soln');
t3 = nexttile;
imagesc(Ubest);title('best soln');
t4 = nexttile;
imagesc(U2);title('original BM3DDEB');
linkaxes([t1 t2 t3 t4]);

figure(688);
imagesc(lambdas,Cs,ee);
title('error matrix');
colormap(jet);colorbar;
xlabel('lambda value');
ylabel('threshold constant');
set(gca,'Xtick',lambdas,'Ytick',Cs);
imcontrast