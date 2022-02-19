N = 64;
d = 8;
xRi = [32 45];

I = zeros(N);
I(xRi(1),xRi(2)) = 1;
% I(1:4,1:4) = 1;

g = zeros(N,N);
g(1:d,1:d) = 1;
blockNorms = real(circshift(ifft2(fft2(g).*fft2(I.^2)),[-d+1,-d+1]));

xR = I(xRi(1):xRi(1)+d-1,xRi(2):xRi(2)+d-1);


FI = fft2(I);
IP = real(ifft2(conj(fft2(xR,N,N)).*FI));
blockDistances = (blockNorms(xRi(1),xRi(2)) + blockNorms - 2*IP);


figure(123);
subplot(2,2,1);imagesc(blockNorms);
subplot(2,2,2);imagesc(IP);
subplot(2,2,3);imagesc(blockDistances);
