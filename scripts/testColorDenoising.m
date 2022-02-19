% testing the color-BM3D to denoise RGB images. It does not seem to improve
% much (if at all) over just denoising each channel one-by-one. The only
% advantage seems that you block match just once instead of 3 times. 


clear;
SNR = 15*(rand(3,1) + .1);
path = 'C:\Users\toby.sanders\Dropbox\archives\data\testImages\';
I = im2double(imread([path,'lena.png']));
b = I;
sigma = zeros(1,1,3);
for i = 1:3
    [b(:,:,i),sigma(i)] = add_Wnoise(I(:,:,i),SNR(i));
end

y1 = CBM3D(b,sigma);
y2 = BM3D(b,sigma);
y3 = zeros(size(b,1),size(b,2),3);
for i = 1:3
    y3(:,:,i) = BM3D(b(:,:,i),sigma(i));
end

%%
fprintf('CBM3D error: %g\n',myrel(I,y1));
fprintf('BM3D error, each channel processed separately: %g\n',myrel(I,y3));
fprintf('alt. M3D error: %g\n',myrel(I,y2));


figure(177);tiledlayout(2,2,'tilespacing','none');
t1 = nexttile;imagesc(b);title('noisy');
t2 = nexttile;imagesc(y1);title('CBM3D');
t3 = nexttile;imagesc(y2);
t4 = nexttile;imagesc(y3);title('BM3D, each channel processed separately');