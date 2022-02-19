clear;
d = 256;
sZ = 8;
SNRs = 2:2:16;
minB = [4 8 20 50 100];
k = 2;
opts.blockSize = sZ;
opts.numMax = 400;
opts.numMin = 8;
lambda = .0367;


hopts.order = 2;
hopts.levels = 1;
hopts.iter = 75;
hopts.wrap_shrink = false;

rng(2021);

nS = numel(SNRs);
nB = numel(minB);
ee = zeros(nS,nB,3);

% I = imresize(im2double(rgb2gray(imread('peppers.png'))),[d,d]);
I = im2double(imread('cameraman.tif'));
% I = im2double(rgb2gray(imread('IMG_0086.JPG')));
% I = I(1201:1201+d-1,1601:1601+d-1);
% I = phantom(d);
IO = I;

for i = 1:nS
    for j = 1:nB
    
        SNR = SNRs(i);
        opts.numMin = minB(j);
        
        
        [I,sigma] = add_Wnoise(IO,SNR);
        hopts.mu = lambda/sigma^2;
        S = matchBlocks(I,sigma,opts);
        dopts.numMax = opts.numMax;
        [U,out] = denoiseBlocks(I,S,sZ,sigma,dopts);

        recBM = BM3D(I,sigma);
        h = zeros(d,d);h(1) = 1;
        hopts.order = 1;
        hopts.mode = 'deconv';
        recTV = HOTV3D(h,I,[d,d,1],hopts);
        ee(i,j,1) = myrel(U(50:200,50:200),IO(50:200,50:200),1);
        ee(i,j,2) = myrel(recBM(50:200,50:200),IO(50:200,50:200),1);
        ee(i,j,3) = myrel(recTV(50:200,50:200),IO(50:200,50:200),1);
    end
end

%%
figure(323);
subplot(2,2,1);imagesc(minB,SNRs,ee(:,:,1));colormap(jet);
subplot(2,2,2);hold off;
ll = {};
for i = 1:nB
    plot(SNRs,ee(:,i,1));
    hold on;
    ll{i} = num2str(minB(i));
end
legend(ll);