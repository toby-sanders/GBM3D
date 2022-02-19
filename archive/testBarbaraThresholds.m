clear;
d = 128; % image dim.
SNR = 5; % SNR of noisy image
opts.numMax = 16; % max number of blocks to match for each ref. block
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
opts.Wiener = true;
opts.filtType = 'ht';
path = 'C:\Users\toby.sanders\Dropbox\archives\data\testImages\';
rng(2021);

% get image and add noise
% I = imresize(im2double(rgb2gray(imread('peppers.png'))),[d,d]);
% I = im2double(imread('cameraman.tif'));
% I = im2double(imread([path,'barbara.png']));
% I = im2double(rgb2gray(imread('images/IMG_0086.JPG')));
% I = I(1201:1201+d-1,1601:1601+d-1);
% I = phantom(d);
% I = im2double(rgb2gray(imread([path,'lena.png'])));
% I = im2double(rgb2gray(imread([path,'fabio.jpeg'])));
I = im2double(rgb2gray(imread([path,'peppers2.png'])));
% I = im2double(rgb2gray(imread([path,'monarch.png'])));
% I = im2double(rgb2gray(imread([path,'tulips.png'])));
% I = im2double(rgb2gray(imread([path,'cat.png'])));
% I = im2double(rgb2gray(imread([path,'baboon.png'])));
% I = ones(d)*1/2;
% I(3*d/8:5*d/8,3*d/8:5*d/8) = 1;
% I = imread('images/house.tif');
% I = im2double(I(:,:,1));
% I = im2double((imread([path,'surfer.png'])));

IO = I;
[I,sigma] = add_Wnoise(I,SNR);

% our BM3D code
gpuID = gpuDevice(1);
reset(gpuID);


% pad image so its dimensions are divisible by the blocksize
[m,n] = size(I);
M = ceil(m/opts.blockSize)*opts.blockSize;
N = ceil(n/opts.blockSize)*opts.blockSize;

% reflect the image to pad to new dimension
I2 = zeros(M,N);
I2(1:m,1:n) = I;
I2(m+1:M,1:n) = flipud(I(m-(M-m)+1:m,1:n));
I2(1:m,n+1:N) = fliplr(I(1:m,n-(N-n)+1:n));
I2(m+1:M,n+1:N) = fliplr(flipud(I(m-(M-m)+1:m,n-(N-n)+1:n)));
U = zeros(M,N,numel(opts.matchSpins));
W = U;
W1 = U;
U1 = U;

% precompute threshold levels tau and wavelet filters
if size(sigma,1)==1 
    tau = sigma*getTau(opts.levels,opts.tauMode);
else
    tau = getColoredTau(opts.levels,opts.order,sigma);
    opts.Wiener = false;
end
if strcmp(opts.wname,'bior')
    [Psi,Psi2] = getWaveFilters3Dbior(opts.wnamez,opts.order,opts.orderz);
else
    Psi = getWaveFilters3D(opts.wname,opts.wnamez,opts.order,opts.orderz);
    Psi2 = Psi;
end



[S,Volume] = matchBlocksCPU(I2,sigma,opts);

nTau = 6;
tauAll = linspace(2.3,3.2,nTau);
ee = zeros(nTau,nTau,nTau);
for i1 = 1:nTau
    for i2 = 1:nTau
        % gpuID = gpuDevice(1);
        % reset(gpuID);
        for i3 = 1:nTau
            tau(1,:) = tauAll(i1);
            tau(2,:) = tauAll(i2);
            tau(3,:) = tauAll(i3);
            tau = tau*sigma
            [U,W] = denoiseAndAggregateGPU(Volume,S,tau,Psi,Psi2,opts);
            U = U./W;
            ee(i1,i2,i3) = myPSNR(IO,U,1)
        end
    end
end


