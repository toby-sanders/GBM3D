clear;
d = 512; % image dim.
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

P = phantom(d);
[P,sigma] = add_Wnoise(P,SNR);

[S1,V1] = matchBlocksCPUnew(P,sigma,opts);

tic;
[S2,V2] = matchBlocksCPU(P,sigma,opts);
toc;