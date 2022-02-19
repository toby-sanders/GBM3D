clear;
sigma = .05;
scl = 4;


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
opts.filtType = 'ht';

T = im2double(rgb2gray(imread('images/SARim.jpg')));

rec = LTBM3D(T,sigma,opts);



%%
dbr = 60;
figure(111);tiledlayout(1,2,'tilespacing','none');colormap(gray)
t1 = nexttile;imagesc(T);
t2 = nexttile;imagesc(rec);
linkaxes([t1 t2]);

figure(112);displaySAR(rec);