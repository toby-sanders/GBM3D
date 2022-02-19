% testing method for determining compression/blocking artifacts blindly

clear;
d = 128; % image dim.
sZ = 32; % size of blocks
SNR = 10; % SNR of noisy image
opts.blockSize = sZ;
opts.numMax = 16; % max number of blocks to match for each ref. block
opts.numMin = 16; % min number of blocks to match for each ref. block
opts.wname = 'bior';
opts.wnamez = 'db';
opts.order = 13;
opts.orderz = 1;
opts.levels = 4;
opts.cycleSpin = 2;
opts.matchSpins =  [0 5 10 20];
opts.wiener = false;
opts.tauMode = 2;
rng(2021);



sL = 128; % length of 1D FT sequences
chooseCase = 5; % choose which data set
Vpath = 'C:\Users\toby.sanders\Dropbox\archives\data\NGA_ResultsAndData\mainTestMovies';
myDir = pwd;

% read in data
if chooseCase==1
    
    V = VideoReader([Vpath,filesep,'ManOnRoof.mpg']);
    slice = 1200;
    yy = 201:800; % image cropping dimensions
    xx = 251:850;
elseif chooseCase==2
    Vpath = 'C:\Users\toby.sanders\Dropbox\archives\data\NGA_ResultsAndData\Test Movies 1-17-20';
    % cd('../Test Movies 1-17-20');
    V = VideoReader([Vpath,filesep,'shippingContainers.mpg']);
    xx = 50:700;
    yy = 50:440;
    slice = 100;
elseif chooseCase==3
    Vpath = 'C:\Users\toby.sanders\Dropbox\archives\data\NGA_ResultsAndData\Test Movies 1-17-20';
     % cd('../Test Movies 1-17-20');
     V = VideoReader([Vpath,filesep,'soccerField.mpg']);
     xx = 50:680;
     yy = 40:440;
     slice = 700+16-20;
elseif chooseCase==4 % simulated case
    Vpath = 'C:\Users\toby.sanders\Dropbox\archives\data\NGA_ResultsAndData\TestMovies-3-31-20';
    V = VideoReader([Vpath,filesep,'trucksZoom.mpg']);
   slice =601;
   xx = 200:1600;
   yy = 201:850;
elseif chooseCase==5 % simulated case
    Vpath = 'C:\Users\toby.sanders\Dropbox\archives\data\NGA_ResultsAndData';
    V = VideoReader([Vpath,filesep,'NighttimeSoccerField.mpg']);
   slice =98;
   xx =1:1900;
   yy = 1:1080;
end



f1 = rgb2gray(readFrame(V));
Vid = zeros(V.Height,V.Width,slice,class(f1));

Vid(:,:,1) = f1;
for i = 2:slice+5
    Vid(:,:,i) = rgb2gray(readFrame(V));
end

I = im2double(Vid(yy,xx,slice));


sigma = 9/255;
opts.matchSpins =  [0,4,16];
U = LTBM3D(I,sigma,opts);
U2 = BM3D(I,sigma);

alpha = 2;
IS1 = min(max(LaplacianSharpen(U,alpha),0),1);
    

% BM3D colored denoising of sharpened image
[m,n] = size(IS1);
T = [0 -1 0; -1 4 -1; 0 -1 0];
Lfilt = 1 + alpha*abs(fft2(T,m,n));
Lsigma = 1;
LsigmaPSD = (Lsigma)^2*abs(Lfilt).^2;
IS = BM3D(IS1,LsigmaPSD);


%%


figure(103);tiledlayout(2,3,'tilespacing','none');colormap(gray);
s1 = nexttile;imagesc(I);title('original frame');
set(gca,'fontsize',16,'fontweight','bold','Ytick',[],'Xtick',[]);
s2 = nexttile;imagesc(U);title('denoised frame');
set(gca,'fontsize',16,'fontweight','bold','Ytick',[],'Xtick',[]);
s3 = nexttile;imagesc(IS,[min(U(:)),max(U(:))]);title('denoised and sharpened frame');
set(gca,'fontsize',16,'fontweight','bold','Ytick',[],'Xtick',[]);
s4 = nexttile;imagesc(I);title('original frame');
set(gca,'fontsize',16,'fontweight','bold','Ytick',[],'Xtick',[]);
s5 = nexttile;imagesc(U);title('denoised frame');
set(gca,'fontsize',16,'fontweight','bold','Ytick',[],'Xtick',[]);
s6 = nexttile;imagesc(IS,[min(U(:)),max(U(:))]);title('denoised and sharpened frame');
set(gca,'fontsize',16,'fontweight','bold','Ytick',[],'Xtick',[]);
linkaxes([s1 s2 s3]);
linkaxes([s4 s5 s6]);


figure(104);tiledlayout(1,2,'tilespacing','none');colormap(gray);
t1 = nexttile;imagesc(U,[min(U(:)),max(U(:))]);title('G-BM3D');
t2 = nexttile;imagesc(U2,[min(U(:)),max(U(:))]);title('BM3D');
linkaxes([t1 t2]);