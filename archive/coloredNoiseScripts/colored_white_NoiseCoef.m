sigma = 0.5;
omega = .0001;
lambda = 1;
d = 128;
wname = 'db';
order = 1;
levels = 3;

e = randn(d)*sigma;
Fe = fft2(e)/d;



Psi = getWaveFilters2D(wname,order);
C = myWavDec2DFFT(Psi,levels,e);

figure(654);
subplot(2,2,1);imagesc(abs(e));colorbar;
subplot(2,2,2);imagesc(abs(Fe));colorbar;
% subplot(2,2,3);imagesc(recWei);colorbar;
% subplot(2,2,4);imagesc(fftshift(filt));colorbar;

s = zeros(3,3);
figure(655);
for j = 1:3
for i = 1:3
    subplot(3,3,(j-1)*3+i);
    imagesc(C{i,j});colorbar;
    s(i,j) = mean(col(abs(C{i,j})));
end
end