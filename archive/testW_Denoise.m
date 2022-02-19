d = 256;
SNR = 5;
levels = 1;

P = phantom(d);
% P = zeros(d,d);
% P(d/2,d/2) = 1;
P = add_Wnoise(P,SNR);
rec = my_wavelet_denoise_3D('db1',levels,P);

figure(234);
subplot(2,2,1);imagesc(P);colorbar;
subplot(2,2,2);imagesc(rec);colorbar;