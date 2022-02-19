levels = 3;
k = 2;
p = 512;
q = 1;
r = 1;
SNR = 3;

l = 0:levels-1;
l = 2.^l;
VX = zeros(p,q,r,levels); VY = VX; VZ = VX;
for ii = 1:levels
    vx = (-1)^k*[0,1/l(ii)*((exp(-1i*2*pi*(1:q-1)*l(ii)/q)-1).^(k+1))./(exp(-1i*2*pi*(1:q-1)/q)-1)];
    vy = (-1)^k*[0,1/l(ii)*((exp(-1i*2*pi*(1:p-1)*l(ii)/p)-1).^(k+1))./(exp(-1i*2*pi*(1:p-1)/p)-1)];
    vz = (-1)^k*[0,1/l(ii)*((exp(-1i*2*pi*(1:r-1)*l(ii)/r)-1).^(k+1))./(exp(-1i*2*pi*(1:r-1)/r)-1)];
    [VX(:,:,:,ii),VY(:,:,:,ii),VZ(:,:,:,ii)] = meshgrid(vx,vy,vz);
end


figure(1);hold off;
F = ones(p,1);
for i = 1:levels
    VY(:,:,:,i) = abs(VY(:,:,:,i))/max(abs(col(VY(:,:,:,i))))/2;
    figure(1);
    plot(abs(VY(:,1,1,i)));
    hold on;
    F = F- VY(:,1,1,i);
end
plot(F);

x = zeros(p,1);
x(p/4:3*p/4) = 1;
x = add_Wnoise(x,SNR);
Fx = fft(x);
xD = zeros(p,levels+1);
xD(:,end) = ifft(Fx.*F);
for i = 1:levels
    xD(:,i) = ifft(VY(:,1,1,i).*Fx);
end

figure(2);hold off;
for i = 1:levels+1
    figure(2);
    plot(real(xD(:,i)));hold on;
end
