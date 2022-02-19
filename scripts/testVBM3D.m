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
opts.matchSpins =  [0 2 4 6 8 10 12 14];
opts.wiener = false;
opts.tauMode = 2;
rng(2021);



sL = 128; % length of 1D FT sequences
chooseCase = 3; % choose which data set
Vpath = 'C:\Users\toby.sanders\Dropbox\archives\data\NGA_ResultsAndData\mainTestMovies';
myDir = pwd;

% read in data
if chooseCase==1
    
    V = VideoReader([Vpath,filesep,'ManOnRoof.mpg']);
    slice = 1000;
    yy = 201:800; % image cropping dimensions
    xx = 251:850;
elseif chooseCase==2
    Vpath = 'C:\Users\toby.sanders\Dropbox\archives\data\NGA_ResultsAndData\Test Movies 1-17-20';
    cd('../Test Movies 1-17-20');
    V = VideoReader('shippingContainers.mpg');
    xx = 50:700;
    yy = 50:440;
    slice = 100;
elseif chooseCase==3
    Vpath = 'C:\Users\toby.sanders\Dropbox\archives\data\NGA_ResultsAndData\Test Movies 1-17-20';
     % cd('../Test Movies 1-17-20');
     V = VideoReader([Vpath,filesep,'soccerField.mpg']);
     xx = 50:680;
     yy = 40:440;
     slice = 6;
elseif chooseCase==4 % simulated case
    SNR = 5;
    I0 = im2double(imread('cameraman.tif'));
    [I,sigma] = add_Wnoise(I0,SNR);
    
    % perform hard thresholding of wavelet coefficients to generate blocking artifacts
    tau = 4*sigma;
    % tau = 0;
    % cd '/Users/tobysanders/Dropbox/TobySharedMATLAB/TobyShared/solvers/DenoisingEngines/LTBM3D'
    Psi = getWaveFilters2DCPU('db',1);
    levels = 4;
    C = myWavDec2DFFT(Psi,levels,I);
    for i = 1:levels
        for j = 1:3
            C{i,j} = C{i,j}.*(abs(C{i,j})>tau);        
        end
    end
    I= myWavRec2D(Psi,levels,C);
    myrel(I,I0)
end

if chooseCase~=4 % in case we're not using the simulated version
    f1 = rgb2gray(readFrame(V));
    Vid = zeros(V.Height,V.Width,slice,class(f1));

    Vid(:,:,1) = f1;
    for i = 2:slice+5
        Vid(:,:,i) = rgb2gray(readFrame(V));
    end

    I = im2double(Vid(yy,xx,slice-5:slice+5));
end

sigma = 7/255;
[U] = VBM3D(I,sigma,opts);
opts.matchSpins =  [0,4,16];
U2 = LTBM3D(I(:,:,6),sigma,opts);

%%
figure(867);colormap(gray);
tiledlayout(2,2,'tilespacing','none');
t1 = nexttile;imagesc(U(:,:,1));title('simple VBM3D');
t2 = nexttile;imagesc(U2);title('BM3D');
t3 = nexttile;imagesc(I(:,:,6));title('noisy');
linkaxes([t1 t2 t3]);