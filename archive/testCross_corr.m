d = 64;
ds = 8;


f = zeros(d,1);
f(1:d/2) = 1;
g = ones(8,1);

fPad = [f;f];
fCg = zeros(d,1);
for i = 1:d
    fCg(i) = sum(fPad(i:i+ds-1).*g);
end


fCg2 = real(ifft(conj(fft(g,d)).*fft(f)));
fCg3 = real(ifft(fft(g,d).*conj(fft(f))));
myrel(fCg,fCg2)
myrel(fCg,fCg3)


figure(121);hold off;
plot(fCg);hold on;
plot(fCg2);
plot(fCg3)
hold off;