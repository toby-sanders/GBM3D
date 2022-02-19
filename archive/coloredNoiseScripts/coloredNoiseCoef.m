clear;
sigma = .1;
omega = 2;
lambda = 10;
d = 256;
wname = 'db';
order = 3;
levels = 3;

e = randn(d)*sigma;
P = phantom(d);
Fe = fft2(e)/d;



[h,hhat] = makeGausPSF(d,omega);
b = real(ifft2(fft2(P).*hhat)) + e;
V = my_Fourier_filters(2,2,d,d,1);
filt = conj(hhat)./(conj(hhat).*hhat + sigma^2*lambda*V);
recWei = real(ifft2(filt.*fft2(e)));
recWei2 = real(ifft2(filt.*fft2(b)));

Psi = getWaveFilters2D(wname,order);
C = myWavDec2DFFT(Psi,levels,recWei);
C2 = myWavDec2DFFT(Psi,levels,recWei2);
figure(654);
subplot(2,2,1);imagesc(abs(e));colorbar;
subplot(2,2,2);imagesc(abs(Fe));colorbar;
subplot(2,2,3);imagesc(recWei);colorbar;
subplot(2,2,4);imagesc(fftshift(filt));colorbar;

s = zeros(3,3);tau = s;
figure(655);
for j = 1:3
for i = 1:3
    figure(655);
    subplot(3,3,(j-1)*3+i);
    imagesc(C{i,j});colorbar;
    s(j,i) = mean(col(abs(C{i,j})));
    tau(j,i) = 3*1.25*s(j,i);
    C2{i,j} = C2{i,j}.*(abs(C2{i,j})>tau(j,i));
    figure(657);
    subplot(3,3,(j-1)*3+i);
    imagesc(C2{i,j});colorbar;
end
end
%%
U = myWavRec2DGPUfast(Psi,levels,C2);
figure(656);
subplot(2,2,1);imagesc(s);colorbar;
subplot(2,2,2);imagesc(recWei2,[0 1]);colorbar;
subplot(2,2,3);imagesc(U,[0 1]);colorbar;
subplot(2,2,4);imagesc(b,[0 1]);colorbar;