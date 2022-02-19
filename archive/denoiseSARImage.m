clear;
sigma0 = .1;
sigma = .001;
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
opts.tauMode = 2;
opts.filtType = 'ht';
opts.Wiener = true;
load('C:\Users\toby.sanders\Dropbox\TobySharedMATLAB\SandBox\radar\gotchaBPimg.mat');
T = abs(T);

rec = LTBM3D(T*scl,sigma,opts)/scl;

logT = (log(T) + 8)/5;
logT = min(max(logT,0),1);
logRec = LTBM3D(logT,sigma0,opts);
opts.filtType = 'ht';
logRec2 = LTBM3D(logT,sigma0,opts);
%%
dbr = 60;
figure(111);tiledlayout(2,2,'tilespacing','none');
t1 = nexttile;image_radar(T,dbr,50);
t2 = nexttile;image_radar(rec,dbr,50);
t3 = nexttile;imagesc(logT);colorbar;
t4 = nexttile;imagesc(logRec);colorbar;
linkaxes([t1 t2]);linkaxes([ t3 t4]);

gam = 0.7;
figure(112);tiledlayout(1,3,'tilespacing','none');
s1 = nexttile;displaySAR(logT,gam);title('SAR image');
set(gca,'Xtick',[],'Ytick',[],'fontsize',18);
s2 = nexttile;displaySAR(logRec,gam);title('denoised SAR image');
set(gca,'Xtick',[],'Ytick',[],'fontsize',18);
s3 = nexttile;displaySAR(logRec2,gam);title('alternative denoising');
set(gca,'Xtick',[],'Ytick',[],'fontsize',18);