B = 16;
B2 = 32;
d = 512;

xR = rand(B,B);
xR2 = rand(B,B,d/B*d/B);

tic;
cnt = 0;
for i = 1:d/B
    for j = 1:d/B
        cnt = cnt+1;
        xR2(:,:,cnt) = xR;
    end
end
t = fft2(xR2,B2,B2);
toc;

tic;
for i = 1:d/B
    for j = 1:d/B
        t = fft2(xR,B2,B2);
    end
end
toc;