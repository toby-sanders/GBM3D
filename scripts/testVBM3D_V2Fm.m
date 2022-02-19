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
opts.levels = 3;
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
     slice = 600;
     sigma = 12/255;
elseif chooseCase==4 % simulated case
    Vpath = 'C:\Users\toby.sanders\Dropbox\archives\data\NGA_ResultsAndData\TestMovies-3-31-20';
    V = VideoReader([Vpath,filesep,'trucksZoom.mpg']);
   slice = 910;
   xx = 200:1600;
   yy = 201:850;
   sigma = 67/255;
elseif chooseCase==5
    xx = 51:1400;
yy = 51:650;
load(['C:\Users\toby.sanders\Dropbox\archives\data\NGA_ResultsAndData\archive\AlignedFrames\',...
    'AlignedFrames-noisy-ManOnRoof-10-23-19-uint8.mat']);
I = im2double(AlignedFrames(yy,xx,40:50));
end

if chooseCase~=5
f1 = rgb2gray(readFrame(V));
Vid = zeros(V.Height,V.Width,slice,class(f1));

Vid(:,:,1) = f1;
for i = 2:slice+5
    Vid(:,:,i) = rgb2gray(readFrame(V));
end
I = im2double(Vid(yy,xx,slice-2:slice+2));
end

ss = (size(I,3)+1)/2;
sigma = 17/255;
opts.matchSpins =  [0,4,16];
% U2 = LTBM3D(I(:,:,6),sigma,opts);
U = LTBM3D(I(:,:,ss),sigma,opts);
U2 = V2FBM3D(I,sigma,opts);


%%
figure(243);colormap(gray);
tiledlayout(2,2,'tilespacing','none');
t1 = nexttile;imagesc(I(:,:,ss));
t2 = nexttile;imagesc(U);
t3 = nexttile;imagesc(U2);
linkaxes([t1 t2 t3]);