% testing how the best reg. lambda for a Wiener filter correpsonds to the 
% best lambda needed for BM3D. In general, it is clear that the lambda for
% BM3D should be smaller than the ideal for Wiener filter, since BM3D comes
% in behind and cleans up the noise. So to maximize the signal from BM3D
% the Wiener filter should leave "some" noise in.

clear;
SNR = 40;
omega = 2;
order = 2;
rng(2021);

% get image and add noise
path = 'C:\Users\toby.sanders\Dropbox\archives\data\testImages\';
% path = '/home/tobysanders/Dropbox/archives/data/testImages/';
I = im2double(rgb2gray(imread([path,'lena.png'])));
% I = im2double((imread([path,'house.tif'])));I = I(:,:,1);
% I = im2double(rgb2gray(imread([path,'peppers2.png'])));
% I = im2double((imread([path,'confFlag.tif'])));
% I = im2double(rgb2gray(imread([path,'PineyMO_Oimage.tif'])));
% I = im2double(rgb2gray(imread([path,'monarch.png'])));
% I = im2double((imread('cameraman.tif')));
% I = phantom(512);
% I = im2double(rgb2gray(imread([path,'IMG_1502.jpeg'])));


[d1,d2] = size(I);
[h,hhat] = makeGausPSF([d1,d2],omega,omega,0,1);
hhat2 = abs(hhat).^2;
b = real(ifft2(fft2(I).*hhat));
[b,sigma] = add_Wnoise(b,SNR);
Fb = fft2(b);


V = my_Fourier_filters(order,1,d1,d2,1);

opts.sigma = sigma;
opts.V = V;
opts.iter = 250;
opts.tol = 1e-10;
[U,out] = SURE_FP_Lambda(Fb,hhat,opts);

NL = 20;
lambda = real(out.lambdas(end));
eLam = log10(real(lambda));
lambdaTest = 10.^linspace(eLam-3,eLam+1,NL);
BMopts.profile = 'default';
ee = zeros(NL,2);
recW = zeros(d1,d2,NL);
recB = zeros(d1,d2,NL);
for i = 1:NL
    lam = lambdaTest(i);
    filt = conj(hhat)./(hhat2 + lam*V);
    sigmaPSD = sigma^2*abs(filt).^2;
    u1 = real(ifft2(Fb.*filt));
    u2 = GBM3D(u1,sigmaPSD,BMopts);
    ee(i,1) = myPSNR(I,u1,1);
    ee(i,2) = myPSNR(I,u2,1);
    recW(:,:,i) = u1;
    recB(:,:,i) = u2;
end

BMopts.alpha = lambda/10;
recBnew = GBM3D_deconv(b,hhat,sigma,BMopts);
BMopts.alpha = 27*sigma^2;
recBnew2 = GBM3D_deconv(b,hhat,sigma,BMopts);

%%
[~,WieBest] = max(ee(:,1));
[~,WieWorst] = min(ee(:,1));
[~,BMBest] = max(ee(:,2));
%%
fprintf('best one step BM3D deconv: %g\n',max(ee(:,2)));
fprintf('BM3D deconv with SURE lambda: %g\n',myPSNR(I,recBnew,1));
fprintf('BM3D deconv with default reg: %g\n',myPSNR(I,recBnew2,1));

figure(791);hold off;
semilogx(lambdaTest,ee(:,1));hold on;
semilogx(lambdaTest,ee(:,2));
plot(ones(2,1)*lambda,[25,35])
hold off;
legend({'Wiener PSNR','BM3D PSNR','SURE lambda'},'location','southeast');
xlabel('lambda');
ylabel('PSNR');


figure(792);tiledlayout(2,3,'tilespacing','none');colormap(gray);
t1 = nexttile;
imagesc(recW(:,:,WieBest),[0 1]);
title('best Wiener');
t2 = nexttile;
imagesc(recW(:,:,WieWorst),[0 1]);
title(['worst Wiener']);
t3 = nexttile;
imagesc(recW(:,:,BMBest),[0 1]);
title('Best Wiener for BM3D');
t4 = nexttile;
imagesc(recB(:,:,WieBest),[0 1]);
t5 = nexttile;
imagesc(recB(:,:,WieWorst),[0 1]);
t6 = nexttile;
imagesc(recB(:,:,BMBest),[0 1]);
linkaxes([t1 t2 t3 t4 t5 t6]);

fprintf('ratio of best lambda for Wiener to best lambda for BM3D: %g\n'...
    ,lambdaTest(WieBest)/lambdaTest(BMBest));

