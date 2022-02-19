clear;
d = 512;
SNR = 5;
levels = 20;
order = 3;
rng(2021);

wname = ['db',num2str(order)];
if strcmp(wname,'db1'), load db1Filters;
elseif strcmp(wname,'db2'), load db2Filters;
elseif strcmp(wname,'db3'), load db3Filters;
end
% x = phantom(d);
x = im2double(imread('cameraman.tif'));
d = size(x,1);
x = add_Wnoise(x,SNR);

U2 = my_wavelet_denoise_3D(wname,levels,x);
hopts.mu = 3;
hopts.order = order;
hopts.mode = 'deconv';
h = zeros(d,d);
h(1) = 1;
rec = HOTV3D(h,x,[d,d,1],hopts);

figure(121);hold off;
subplot(2,2,1);imagesc(x);
% subplot(2,2,2);imagesc(U,'-','linewidth',2);
subplot(2,2,3);imagesc(rec);
subplot(2,2,4);imagesc(real(U2));
hold off;

% myrel(x,U)
myrel(x,U2)

% figure(234);hold off;
% plot(real(loF));hold on;
% plot(real(hiF));