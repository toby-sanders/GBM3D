clear;
d = 2560;
SNR = 5;
levels = 3;
order = 3;
tau = .5;
rng(2021);

wname = ['db',num2str(order)];
if strcmp(wname,'db1'), load db1Filters;
elseif strcmp(wname,'db2'), load db2Filters;
elseif strcmp(wname,'db3'), load db3Filters;
end
loF = fft(Lo_d',d);
hiF = fft(Hi_d',d);

x = phantom(d);
x0 = x;
x = add_Wnoise(x,SNR);

U = my_wavelet_denoise_2D_ortho(wname,levels,x,tau);

% tau = .4;
% mm = 16;
% [c,l] = wavedec(x,levels,wname);
% c(mm+1:end) = c(mm+1:end).*(abs(c(mm+1:end))>tau);
% wrec = waverec(c,l,wname);


hopts.mu = 3;
hopts.order = order;
hopts.mode = 'deconv';
h = zeros(d,1);
h(1) = 1;
% rec = HOTV3D(h,x,[d,d,1],hopts);

figure(121);hold off;
subplot(2,2,1);imagesc(x,[0 1]);
subplot(2,2,2);imagesc(U,[0 1]);
% subplot(2,2,3);imagesc(rec);
hold off;

myrel(U,x)