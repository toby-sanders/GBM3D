clear;
d = 128;
SNR = 15;
levels = 3;
order = 1;
rng(2021);
tau = 0;

wname = ['db',num2str(order)];
if strcmp(wname,'db1'), load db1Filters;
elseif strcmp(wname,'db2'), load db2Filters;
elseif strcmp(wname,'db3'), load db3Filters;
end
loF = fft(Lo_d',d);
hiF = fft(Hi_d',d);

x = zeros(d,1);
x = sin(1.5*pi*linspace(0,1,d)');
x(d/4:3*d/4) = 1;
x = add_Wnoise(x,SNR);

% tmp = x;
% tmp2 = tmp;
% for i = 1:levels
%     tmp = imfilter(tmp,Lo_d');
%     tmp2 = ifft(fft(tmp2).*conj(loF));
%     figure(652);hold off;
%     plot(tmp);hold on;
%     plot(tmp2);
%     pause;
% end



U = my_wavelet_denoise_1D(wname,levels,x,tau);
U2 = my_wavelet_denoise_3D(wname,levels,x);
hopts.mu = 3;
hopts.order = order;
hopts.mode = 'deconv';
h = zeros(d,1);
h(1) = 1;
rec = HOTV3D(h,x,[d,1,1],hopts);

figure(121);hold off;
subplot(2,2,1);plot(x,'linewidth',2);
subplot(2,2,2);plot(U,'-','linewidth',2);
subplot(2,2,3);plot(rec,'linewidth',2);
subplot(2,2,4);plot(real(U2),'linewidth',2);
hold off;

myrel(x,U)
myrel(x,U2)

% figure(234);hold off;
% plot(real(loF));hold on;
% plot(real(hiF));