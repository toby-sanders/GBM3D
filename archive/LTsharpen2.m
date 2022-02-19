function I2 = LTsharpen(I,alpha,sigma)



% BM3D colored denoising of sharpened image
[m,n] = size(I);
T = [0 -1 0; -1 4 -1; 0 -1 0];
lambda = 1e-5;

% T = [-1 -1 -1; -1 8 -1; -1 -1 -1];
T2 = zeros(m,n);
T2(1:3,1:3) = T;
Lfilt = 1 + alpha*fft2(circshift(T2,[-1 -1]),m,n);


a = Lfilt;
b = -1;
c = lambda*Lfilt;

hhat = (-b+sqrt(b^2 - 4.*a.*c))./(2*a);
h = ifft2(hhat);


if numel(sigma)==1
    sigma = ones(m,n)*sigma^2;
end


opts.profile = 'default';
I2 = GBM3D_deconv(I,hhat,sigma,opts);


% ppp = .01;
% figure(111);tiledlayout(2,2,'tilespacing','none');colormap(gray);
% t1 = nexttile;[mm1,mm2] = myImage(I,ppp);
% t2 = nexttile;imagesc(I1,[mm1,mm2]);
% t3 = nexttile;imagesc(I2,[mm1,mm2]);
% t4 = nexttile;imagesc(fftshift(PSD));colorbar;
% linkaxes([t1 t2 t3]);

